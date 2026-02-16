library(nimble)

nimble_model_code <- nimbleCode({

  # Item parameters
  for (j in 1:J) {
    g1[j] ~ dbeta(1, 1)
    g2[j] ~ dbeta(1, 1)
    g3[j] ~ dbeta(1, 1)
    
    s1[j] ~ dbeta(1, 1)
    s2[j] ~ dbeta(1, 1)
    s3[j] ~ dbeta(1, 1)
  }
  
  # Regression for latent transition
  ## Initial latent state
  for (k in 1:K) {
    beta0[k] ~ dnorm(0, 1)
    
    for (p in 1:P) {
      betaZ[k, p] ~ dnorm(0, 1)
    }
    
    ## Transition parameter from non-mastery to mastery
    gamma01[k, 1] ~ dnorm(0, 1)
    #gamma10[k,1]  ~ dnorm(0, 1) # mastery to non-mastery cannot happen
    
    for (p in 1:P) {
      gamma01[k, p + 1] ~ dnorm(0, 1)
    }
  }
  
  # Latent attributes
  for (n in 1:N) {
    for (k in 1:K) {
      logit(p1[n, k]) <- beta0[k] + inprod(betaZ[k, 1:P], Z[n, 1:P])
      alpha1[n, k] ~ dbern(p1[n, k])
    }
    
    ## transition form t1 to t2
    for (k in 1:K) {
      logit(p01_t12[n, k]) <- gamma01[k, 1] + inprod(gamma01[k, 2:(P + 1)], Z[n, 1:P])
      prob_alpha2[n, k] <- alpha1[n, k] + (1 - alpha1[n, k]) * p01_t12[n, k]
      alpha2[n, k] ~ dbern(prob_alpha2[n, k])
    }
    
    ## transition from t2 to t3
    for (k in 1:K) {
      logit(p01_t23[n, k]) <- gamma01[k, 1] + inprod(gamma01[k, 2:(P + 1)], Z[n, 1:P])
      prob_alpha3[n, k] <- alpha2[n, k] + (1 - alpha2[n, k]) * p01_t23[n, k]
      alpha3[n, k] ~ dbern(prob_alpha3[n, k])
    }
  }
  
  # Measurement DINA model
  for (n in 1:N) {
    for (j in 1:J) {
      
      ## time 1
      for (k in 1:K) {
        a_pow1[n, j, k] <- pow(alpha1[n, k], Q1[j, k])
      }
      eta1[n, j] <- prod(a_pow1[n, j, 1:K])
      pY1[n, j]  <- (1 - s1[j]) * eta1[n, j] + g1[j] * (1 - eta1[n, j])
      Y1[n, j]   ~ dbern(pY1[n, j])
      
      ## time 2
      for (k in 1:K) {
        a_pow2[n, j, k] <- pow(alpha2[n, k], Q2[j, k])
      }
      eta2[n, j] <- prod(a_pow2[n, j, 1:K])
      pY2[n, j]  <- (1 - s2[j]) * eta2[n, j] + g2[j] * (1 - eta2[n, j])
      Y2[n, j]   ~ dbern(pY2[n, j])
      
      ## time 3
      for (k in 1:K) {
        a_pow3[n, j, k] <- pow(alpha3[n, k], Q3[j, k])
      }
      eta3[n, j] <- prod(a_pow3[n, j, 1:K])
      pY3[n, j]  <- (1 - s3[j]) * eta3[n, j] + g3[j] * (1 - eta3[n, j])
      Y3[n, j]   ~ dbern(pY3[n, j])
    }
  }
  
  # Posterior probabilities
  for (n in 1:N) {
    for (k in 1:K) {
      prob_attr1[n, k] <- p1[n, k]
      prob_attr2[n, k] <- prob_alpha2[n, k]
      prob_attr3[n, k] <- prob_alpha3[n, k]
    }
  }
})

run_one_simulation_nimble_unknownQ <- function(seed_number = 123,
                                               N = 200,
                                               J = 6,
                                               Q_matrix = Q_matrix,
                                               n.chains = 2,
                                               n.adapt  = 1000,   # nimble这里保留接口但不强制使用
                                               n.burnin = 1000,
                                               n.iter   = 2000,
                                               thin     = 1) {
  
  set.seed(seed_number)
  K <- 3
  P <- 6
  
  # Data Generation (same as yours)
  Z <- MASS::mvrnorm(N, mu = rep(0, P), Sigma = diag(P))
  Z <- scale(Z)
  
  ## True Q-matrices
  Q1 <- matrix(
    c(
      1,0,0,
      0,1,0,
      0,0,1,
      1,1,0,
      1,0,1,
      0,1,1,
      1,0,0,
      0,1,0,
      0,0,1,
      1,0,0,
      0,1,0,
      0,0,1,
      1,0,0,
      0,1,0,
      0,0,1,
      1,0,0,
      0,1,0,
      1,1,0
    ), nrow = 18, byrow = TRUE)
  
  Q2 <- matrix(
    c(
      1,0,0,
      0,1,0,
      0,0,1,
      1,1,0,
      1,0,1,
      0,1,1,
      0,1,0,
      1,0,0,
      0,0,1,
      0,1,0,
      1,0,0,
      0,0,1,
      0,1,0,
      1,0,0,
      0,0,1,
      1,0,0,
      0,1,0,
      1,0,1
    ), nrow = 18, byrow = TRUE)
  
  Q3 <- matrix(
    c(
      1,0,0,
      0,1,0,
      0,0,1,
      1,1,0,
      1,0,1,
      0,1,1,
      1,0,0,
      0,1,0,
      0,0,1,
      1,0,0,
      0,1,0,
      0,0,1,
      1,0,0,
      0,1,0,
      0,0,1,
      0,1,1,
      1,0,0,
      1,0,0
    ), nrow = 18, byrow = TRUE)
  
  Q1_true <- Q1
  Q2_true <- Q2
  Q3_true <- Q3

  
  ## True item parameters
  gs_true_t1 <- matrix(runif(J * 2, 0.05, 0.20), ncol = 2)
  gs_true_t2 <- matrix(runif(J * 2, 0.05, 0.20), ncol = 2)
  gs_true_t3 <- matrix(runif(J * 2, 0.05, 0.20), ncol = 2)
  
  ## True initial attribute states
  beta0 <- c(-1.12, -1.193, -1.00)
  betaZ <- matrix(
    c(
      -0.576,  0.092,  1.282, -0.03,  -0.196,  0.06,
      0.125, -0.302, -0.03,   0.315, -0.567,  0.261,
      0.20,   0.10,   0.30,   0.00,  -0.10,   0.05
    ),
    nrow = K, byrow = TRUE
  )
  
  ## True transition parameters
  gamma01 <- matrix(
    c(
      -0.021, -0.681, -0.192,  0.083,  0.379, -0.104, -0.075,
      0.024, -0.054,  0.109,  0.191,  0.059,  0.125, -0.456,
      0.10,  -0.20,   0.05,   0.15,  -0.10,   0.05,  -0.30
    ),
    nrow = K, byrow = TRUE
  )
  
  invlogit <- function(x) 1 / (1 + exp(-x))
  
  ## Attribute profiles
  p_init <- t(apply(Z, 1, function(z) invlogit(beta0 + betaZ %*% z)))
  alpha1 <- matrix(rbinom(N * K, 1, as.vector(p_init)), N, K)
  
  alpha2 <- matrix(0, N, K)
  alpha3 <- matrix(0, N, K)
  
  for (i in 1:N) {
    z_ext <- c(1, Z[i, ])
    for (k in 1:K) {
      p01_t12 <- invlogit(sum(gamma01[k, ] * z_ext))
      alpha2[i, k] <- if (alpha1[i, k] == 1) 1 else rbinom(1, 1, p01_t12)
      
      p01_t23 <- invlogit(sum(gamma01[k, ] * z_ext))
      alpha3[i, k] <- if (alpha2[i, k] == 1) 1 else rbinom(1, 1, p01_t23)
    }
  }
  
  gen_Y <- function(alpha, Q, gs) {
    outer(seq_len(N), seq_len(nrow(Q)), Vectorize(function(i, j) {
      eta <- prod(alpha[i, ] ^ Q[j, ])
      p   <- (1 - gs[j, 2]) * eta + gs[j, 1] * (1 - eta)
      rbinom(1, 1, p)
    }))
  }
  
  Y1 <- gen_Y(alpha1, Q1_true, gs_true_t1)
  Y2 <- gen_Y(alpha2, Q2_true, gs_true_t2)
  Y3 <- gen_Y(alpha3, Q3_true, gs_true_t3)
  
  # nimble inputs (same info as jags_data but split)
  constants <- list(N = N, J = J, K = K, P = P,Q1 = Q1_true, Q2 = Q2_true, Q3 = Q3_true)
  
  data_list <- list(Z = Z, Y1 = Y1, Y2 = Y2, Y3 = Y3)
  
  # params monitor (same as yours)
  params_monitor <- c("g1","s1","g2","s2","g3","s3",
                      "beta0","betaZ","gamma01",
                      "alpha1","alpha2","alpha3")
  
  # Initial values (same logic)
  make_inits <- function() {
   
    list(
      g1 = runif(J, 0.05, 0.1), s1 = runif(J, 0.05, 0.1),
      g2 = runif(J, 0.05, 0.1), s2 = runif(J, 0.05, 0.1),
      g3 = runif(J, 0.05, 0.1), s3 = runif(J, 0.05, 0.1),
      
      beta0 = rnorm(K, -1.1, 0.05),
      betaZ = matrix(rnorm(K * P, -0.2, 0.2), K, P),
      gamma01 = matrix(rnorm(K * (P + 1), 0, 0.1), K, (P + 1))
    )
  }
  
  inits_list <- replicate(n.chains, make_inits(), simplify = FALSE)
  
  # build nimble model
  Rmodel <- nimbleModel(code = nimble_model_code,
                        constants = constants,
                        data = data_list,
                        inits = inits_list[[1]])  # 先给一个，用于建模结构
  
  Cmodel <- compileNimble(Rmodel)
  
  # MCMC configuration (不改你逻辑：默认采样器，让nimble自己分配)
  conf <- configureMCMC(Rmodel, monitors = params_monitor, thin = thin)
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  # runMCMC: multiple chains
  samples <- runMCMC(Cmcmc,
                     niter = n.burnin + n.iter,
                     nburnin = n.burnin,
                     nchains = n.chains,
                     inits = inits_list,
                     samplesAsCodaMCMC = TRUE,   # 让它像 rjags 一样是 coda mcmc.list
                     thin = thin,
                     setSeed = seed_number + 0:(n.chains - 1))
  
  smat <- as.matrix(samples)
   
  list(
    seed = seed_number,
    samples = samples,
    
    true = list(
      Q1_true = Q1_true, Q2_true = Q2_true, Q3_true = Q3_true,
      gs_true_t1 = gs_true_t1, gs_true_t2 = gs_true_t2, gs_true_t3 = gs_true_t3,
      beta0 = beta0, betaZ = betaZ, gamma01 = gamma01,
      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3,
      N = N, J = J, K = K
    )
  )
}

out <- run_one_simulation_nimble_unknownQ(
  seed_number = 1,
  N = 200,
  J = 18,
  n.chains = 2,
  n.burnin = 2000,
  n.iter   = 2000,
  thin     = 5
)

str(out, max.level = 1)
out$success
out$prop_keep_all
nrow(out$smat_keep)



### 6. Convergence diagnostics -------------------------------
diag_mcmc <- function(samples, rhat_thr = 1.10) {
  
  ## 1) Rhat (Gelman-Rubin)
  psrf <- tryCatch(
    gelman.diag(samples, multivariate = FALSE)$psrf[, 1],
    error = function(e) NA
  )
  
  ## 2) exclude individual level alpha
  keep_names <- names(psrf)[!grepl("^alpha[123]\\[", names(psrf))]
  psrf_s <- psrf[keep_names]
  
  ## 3) ESS
  ess <- tryCatch(
    effectiveSize(samples),
    error = function(e) NA
  )
  ess_s <- ess[intersect(names(ess), keep_names)]
  
  ## 4) Summary
  diag_df <- data.frame(
    max_rhat    = max(psrf_s, na.rm = TRUE),
    median_rhat = median(psrf_s, na.rm = TRUE),
    num_rhat_bad = sum(psrf_s > rhat_thr, na.rm = TRUE),
    
    min_ess     = min(ess_s, na.rm = TRUE),
    median_ess  = median(ess_s, na.rm = TRUE),
    mean_ess    = mean(ess_s, na.rm = TRUE)
  )
  
  ## 5) Debug: parameters with high Rhat
  bad_params <- data.frame(
    parameter = names(psrf_s)[psrf_s > rhat_thr],
    rhat      = as.numeric(psrf_s[psrf_s > rhat_thr])
  )
  
  list(
    summary  = diag_df,
    bad_rhat = bad_params
  )
}

diag_results <- diag_mcmc(out$samples)
diag_results$summary
diag_results$bad_rhat

library(coda)

# tract plots
par_trace <- c(
  "g1[1]","s1[1]","g2[1]","s2[1]","g3[1]","s3[1]")

traceplot(out$samples[, par_trace], main = "Trace plots: structural parameters")

### 7. Evaluation (based on smat_keep only) -------------------------------

## ---- 7a. helpers -------------------------------------------------------
get_draws <- function(smat, pattern) {
  cols <- grep(pattern, colnames(smat))
  if (length(cols) == 0) stop("No columns matched pattern: ", pattern)
  smat[, cols, drop = FALSE]
}

get_posterior_mean_vec <- function(smat, pattern) {
  draws <- get_draws(smat, pattern)
  colMeans(draws)
}

get_posterior_mean_mat <- function(smat, pattern, nrow, ncol, byrow = FALSE) {
  v <- get_posterior_mean_vec(smat, pattern)
  matrix(v, nrow = nrow, ncol = ncol, byrow = byrow)
}

evaluate_continuous_matrix <- function(est, true) {
  diff <- est - true
  rmse <- sqrt(mean(diff^2))
  bias <- mean(diff)
  mab  <- mean(abs(diff))  # mean absolute bias
  list(rmse = rmse, bias = bias, mab = mab)
}
evaluate_binary_matrix <- function(est, true) {
  diff <- est - true
  accuracy <- mean(round(est) == true)
  rmse     <- sqrt(mean(diff^2))
  mab      <- mean(abs(diff))
  list(accuracy = accuracy, rmse = rmse, mab = mab)
}

get_ci_mat <- function(smat, pattern) {
  draws <- get_draws(smat, pattern)
  rbind(
    lower = apply(draws, 2, quantile, 0.025),
    mean  = colMeans(draws),
    upper = apply(draws, 2, quantile, 0.975)
  )
}

compute_PAR <- function(pred, true) mean(apply(pred == true, 1, all))
compute_AAR <- function(pred, true) colMeans(pred == true)

## PIP for alpha (posteriro mastery probabilities)
## Each user have a PIP
alpha_diag_by_user <- function(alpha_true, postprob){
  data.frame(
    mean_pip_on_ones  = rowMeans(postprob*(alpha_true==1), na.rm=TRUE),
    mean_pip_on_zeros = rowMeans(postprob*(alpha_true==0), na.rm=TRUE),
    prop_pip_ambig    = rowMeans(abs(postprob-0.5)<0.05) # hard to decide, range is 0.45-0.55
  )
}

## Evaluation for continuous vector/matrix parameters
eval_cont_vec <- function(est, true) evaluate_continuous_matrix(as.numeric(est), as.numeric(true))
eval_cont_mat <- function(est, true) evaluate_continuous_matrix(as.numeric(est), as.numeric(true))


## ---- 7b. main evaluation function --------------------------------------
evaluate_one_run <- function(smat_keep,
                             J, K, N,
                             Q1_true, Q2_true, Q3_true,
                             gs_true_t1, gs_true_t2, gs_true_t3,
                             beta0, betaZ, gamma01,
                             alpha1, alpha2, alpha3,
                             attr_thr  = 0.5) {
  
 ## 2) Item parameters (posterior means)
  g1_est <- get_posterior_mean_vec(smat_keep, "^g1\\[")
  s1_est <- get_posterior_mean_vec(smat_keep, "^s1\\[")
  g2_est <- get_posterior_mean_vec(smat_keep, "^g2\\[")
  s2_est <- get_posterior_mean_vec(smat_keep, "^s2\\[")
  g3_est <- get_posterior_mean_vec(smat_keep, "^g3\\[")
  s3_est <- get_posterior_mean_vec(smat_keep, "^s3\\[")
  
  ## 3) Regression parameters
  beta0_est <- get_posterior_mean_vec(smat_keep, "^beta0\\[")          # length K
  betaZ_est <- get_posterior_mean_mat(smat_keep, "^betaZ\\[", K, 6)    # K x P (P=6)
  gamma01_est <- get_posterior_mean_mat(smat_keep, "^gamma01\\[", K, 7)# K x (P+1)
  
  ## 4) Latent attributes: posterior prob + hard classification
  a1_postprob <- get_posterior_mean_mat(smat_keep, "^alpha1\\[", N, K)
  a2_postprob <- get_posterior_mean_mat(smat_keep, "^alpha2\\[", N, K)
  a3_postprob <- get_posterior_mean_mat(smat_keep, "^alpha3\\[", N, K)
  a1_hat <- 1L * (a1_postprob >= attr_thr)
  a2_hat <- 1L * (a2_postprob >= attr_thr)
  a3_hat <- 1L * (a3_postprob >= attr_thr)
  
  ## 5) Recovery metrics
  res_g1 <- eval_cont_vec(g1_est, gs_true_t1[, 1])
  res_s1 <- eval_cont_vec(s1_est, gs_true_t1[, 2])
  res_g2 <- eval_cont_vec(g2_est, gs_true_t2[, 1])
  res_s2 <- eval_cont_vec(s2_est, gs_true_t2[, 2])
  res_g3 <- eval_cont_vec(g3_est, gs_true_t3[, 1])
  res_s3 <- eval_cont_vec(s3_est, gs_true_t3[, 2])
  
  ## 5c) regression parameters
  res_beta0   <- eval_cont_vec(beta0_est, beta0)
  res_betaZ   <- eval_cont_mat(betaZ_est, betaZ)
  res_gamma01 <- eval_cont_mat(gamma01_est, gamma01)
  
  ## 5d) attribute classification accuracy
  PAR_t1 <- compute_PAR(a1_hat, alpha1)
  PAR_t2 <- compute_PAR(a2_hat, alpha2)
  PAR_t3 <- compute_PAR(a3_hat, alpha3)
  
  AAR_t1 <- compute_AAR(a1_hat, alpha1)
  AAR_t2 <- compute_AAR(a2_hat, alpha2)
  AAR_t3 <- compute_AAR(a3_hat, alpha3)
  
  diag1 <- alpha_diag_by_user(alpha1, a1_postprob)
  diag2 <- alpha_diag_by_user(alpha2, a2_postprob)
  diag3 <- alpha_diag_by_user(alpha3, a3_postprob)
  # diag1[17,]
  #a1_postprob[17,]
  
  ## 6) Uncertainty: 95% CIs for key continuous parameters
  ci <- list(
    g1 = get_ci_mat(smat_keep, "^g1\\["),
    s1 = get_ci_mat(smat_keep, "^s1\\["),
    g2 = get_ci_mat(smat_keep, "^g2\\["),
    s2 = get_ci_mat(smat_keep, "^s2\\["),
    g3 = get_ci_mat(smat_keep, "^g3\\["),
    s3 = get_ci_mat(smat_keep, "^s3\\["),
    beta0  = get_ci_mat(smat_keep, "^beta0\\["),
    betaZ  = get_ci_mat(smat_keep, "^betaZ\\["),
    gamma01 = get_ci_mat(smat_keep, "^gamma01\\[")
  )
  
  list(
    ## 1) estimates for key parameters
    est = list(
        ## item parameters
      g1 = g1_est, s1 = s1_est,
      g2 = g2_est, s2 = s2_est,
      g3 = g3_est, s3 = s3_est,
      
      ## regressions
      beta0 = beta0_est,
      betaZ = betaZ_est,
      gamma01 = gamma01_est,
      
      ## latent attributes
      a1_postprob = a1_postprob, a2_postprob = a2_postprob, a3_postprob = a3_postprob,
      a1_hat = a1_hat, a2_hat = a2_hat, a3_hat = a3_hat
    ),
    
    ## 2) recovery / evaluation metrics
    metrics = list(
         ## item parameters
      g1 = res_g1, s1 = res_s1,
      g2 = res_g2, s2 = res_s2,
      g3 = res_g3, s3 = res_s3,
      
      ## regression parameters
      beta0 = res_beta0,
      betaZ = res_betaZ,
      gamma01 = res_gamma01,
      
      ## attribute classification
      PAR = c(t1 = PAR_t1, t2 = PAR_t2, t3 = PAR_t3),
      AAR = rbind(t1 = AAR_t1, t2 = AAR_t2, t3 = AAR_t3),
      
      ## individual attribute profile diagnostics (select few to check)
      alpha_diag = list(t1 = diag1, t2 = diag2, t3 = diag3)
    ),
    
    ## 3) uncertainty summaries (95% CIs)
    ci = ci
  )
}

smat_keep <- as.matrix(out$samples)
res1 <- evaluate_one_run(smat_keep,
  J = out$true$J, K = out$true$K, N = out$true$N,
  gs_true_t1 = out$true$gs_true_t1, gs_true_t2 = out$true$gs_true_t2, gs_true_t3 = out$true$gs_true_t3,
  beta0 = out$true$beta0, betaZ = out$true$betaZ, gamma01 = out$true$gamma01,
  alpha1 = out$true$alpha1, alpha2 = out$true$alpha2, alpha3 = out$true$alpha3,
  attr_thr  = 0.5
)


res1
