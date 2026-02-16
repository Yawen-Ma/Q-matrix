
library(nimble)

J <- 18
K <- 3
N <- 200

### 1. Nimble model -------------------------------------------------------
# The main nimble model
dina_lta_model <- nimbleCode({
 
   # Q-matrix
   # hierarchical Bernoulli-Beta prior
 #  theta ~ dbeta(8, 12) # sparsity
  theta ~ dbeta(1, 29) # to encourage sparsity
  
  for (j in 1:J) {
    for (k in 1:K) {
      Q1_est[j, k] ~ dbern(theta)
      Q2_est[j, k] ~ dbern(theta)
      Q3_est[j, k] ~ dbern(theta)
    }
  }
  
  # Item parameters
  for (j in 1:J) {
    g1[j] ~ dbeta(1, 1)
    g2[j] ~ dbeta(1, 1)
    g3[j] ~ dbeta(1, 1)
    
    s1[j] ~ dbeta(1, 1)
    s2[j] ~ dbeta(1, 1)
    s3[j] ~ dbeta(1, 1)
  }
  
  # Covariate effects on transition 
  ## Initial latent state
  for (k in 1:K) {
    beta0[k] ~ dnorm(0, 1)
    for (p in 1:P) {
      betaZ[k, p] ~ dnorm(0, 1)
    }
    
    ## Transition parameter from non-mastery to mastery
    gamma01_12[k, 1] ~ dnorm(0, 1) # t1 to t2
    gamma01_23[k, 1] ~ dnorm(0, 1) # t2 to t3
    #gamma10[k,1]  ~ dnorm(0, 1) # mastery to non-mastery cannot happen
    
    for (p in 1:P) {
     gamma01_12[k, p + 1] ~ dnorm(0, 1)
     gamma01_23[k, p + 1] ~ dnorm(0, 1)
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
      logit(p01_t12[n, k]) <- gamma01_12[k, 1] + inprod(gamma01_12[k, 2:(P + 1)], Z[n, 1:P])
      prob_alpha2[n, k] <- alpha1[n, k] + (1 - alpha1[n, k]) * p01_t12[n, k]
      alpha2[n, k] ~ dbern(prob_alpha2[n, k])
    }
    
    ## transition from t2 to t3
    for (k in 1:K) {
      logit(p01_t23[n, k]) <- gamma01_23[k, 1] + inprod(gamma01_23[k, 2:(P + 1)], Z[n, 1:P])
      prob_alpha3[n, k] <- alpha2[n, k] + (1 - alpha2[n, k]) * p01_t23[n, k]
      alpha3[n, k] ~ dbern(prob_alpha3[n, k])
    }
  }
  
  # Measurement DINA model
  for (n in 1:N) {
    for (j in 1:J) {
      ## time 1
      for (k in 1:K) {
        a_pow1[n,j,k] <- 1 - Q1_est[j,k] + Q1_est[j,k] * alpha1[n,k]
      }
      eta1[n, j] <- prod(a_pow1[n, j, 1:K])
      pY1[n, j]  <- (1 - s1[j]) * eta1[n, j] + g1[j] * (1 - eta1[n, j])
      Y1[n, j]   ~ dbern(pY1[n, j])
      
      ## time 2
      for (k in 1:K) {
        a_pow2[n,j,k] <- 1 - Q2_est[j,k] + Q2_est[j,k] * alpha2[n,k]
      }
      eta2[n, j] <- prod(a_pow2[n, j, 1:K])
      pY2[n, j]  <- (1 - s2[j]) * eta2[n, j] + g2[j] * (1 - eta2[n, j])
      Y2[n, j]   ~ dbern(pY2[n, j])
      
      ## time 3
      for (k in 1:K) {
        a_pow3[n,j,k] <- 1 - Q3_est[j,k] + Q3_est[j,k] * alpha3[n,k]
      }
      eta3[n, j] <- prod(a_pow3[n, j, 1:K])
      pY3[n, j]  <- (1 - s3[j]) * eta3[n, j] + g3[j] * (1 - eta3[n, j])
      Y3[n, j]   ~ dbern(pY3[n, j])
    }
  }
    # Posterior probabilities
 # for (n in 1:N) {
 #   for (k in 1:K) {
  #    prob_attr1[n, k] <- p1[n, k]
  #    prob_attr2[n, k] <- prob_alpha2[n, k]
 #     prob_attr3[n, k] <- prob_alpha3[n, k]
 #   }
#  }
})

library(nimble)
sampler_Q_rowwise <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    
    j     <- control$j
    J     <- control$J
    K     <- control$K
    qname <- control$qname # "Q1_est" or "Q2_est" or "Q3_est"
    
    # candidates: all proposed Q rows (non-zero)
    candidate <- rbind(
      c(0,0,1),
      c(0,1,0),
      c(0,1,1),
      c(1,0,0),
      c(1,0,1),
      c(1,1,0),
      c(1,1,1)
    )
    
    nCandidate <- 7L
    # dependencies (likelihood + priors affected by this row) Y likelihood, alpha likelihood, Q prior, theta
    deps <- model$getDependencies(target, self = FALSE)
  },
  
run = function (){
  Qcurrent <- model[[qname]][1:J, 1:K]
  
  oldRow <- Qcurrent[j, 1:K]
  
  logprob <- nimNumeric(nCandidate, init = TRUE, value = -Inf)
  
  for (m in 1:nCandidate){
    ## always reset model row to oldRow before evaluating candidate m
    for (kk in 1:K) model[[qname]][j, kk] <<- oldRow[kk]
  
    Qstar <- Qcurrent
    
    # candidate row
    for (kk in 1:K) Qstar[j, kk] <- candidate[m, kk]
    
    # check identifiability for the proposed candidate Q row
    
    identifiable <- 1L # nimble is strict to integer, rather than identifiable <- 1 typeof(identifiable) which is double

      # 0) no 0 row
      # if (any(rowSums(Q) == 0)) return(FALSE)
      # 1) each column at least 3 ones
      # if (any(colSums(Q) < 3)) return(FALSE)
      # 2) Q contains an identity submatrix I_K (i.e., at least one pure item e_k for each k)
      #  ek_rows_by_k <- lapply(1:K, function(k) which(rowSums(Q) == 1 & Q[, k] == 1))
      #  if (any(vapply(ek_rows_by_k, length, integer(1)) == 0)) return(FALSE)
      # pick one pure item per attribute
      #   I_idx <- vapply(1:K, function(k) ek_rows_by_k[[k]][1], integer(1))
      # 3) excluding I_K, remaining (J-K) x K has K mutually distinct column vectors
      #   if (self$J <= K) return(FALSE)
      # Q_rem <- Q[-I_idx, , drop = FALSE]
      # col_patterns <- apply(Q_rem, 2, paste0, collapse = "")
      #   if (length(unique(col_patterns)) < K) return(FALSE)
      # (0) no zero row
      for (jj in 1:J) {
        rs <- 0L
        for (kk in 1:K) rs <- rs + Qstar[jj, kk]
        if (rs == 0L) identifiable <- 0L
      }
      if (identifiable == 1L) {
      # (1) each column has at least 3 ones
      for (kk in 1:K) {
        cs <- 0L
        for (jj in 1:J) cs <- cs + Qstar[jj, kk]
        if (cs < 3L) identifiable <- 0L
      }
      }
      if (identifiable == 1L) {
      # (2) at least one pure item for each skill
      # pure item: row sum == 1 and Q[j,k] == 1
      for (kk in 1:K) {
        has_pure <- 0L
        for (jj in 1:J) {
          rs <- 0L
          for (tt in 1:K) rs <- rs + Qstar[jj, tt]
          if (rs == 1L) {
            if (Qstar[jj, kk] == 1L) {
              has_pure <- 1L
            }
          }
        }
        if (has_pure == 0L) identifiable <- 0L
       }
      }
      if(identifiable == 1L){
        
      ## write candidate row into model
      for (kk in 1:K) model[[qname]][j,kk] <<- candidate[m,kk]
      
      ## compute log prob for this candidate: likelihood x priors
      logprob[m] <- model$calculate(deps)
      }
  }
  
  mx <- max(logprob)
  w <- exp(logprob - mx)
  s <- sum(w)
  
  bad <- 0L
  if (s <= 0) bad <- 1L
  if (!(s == s)) bad <- 1L   # NaN check
  
  if (bad == 1L) {
    for (kk in 1:K) model[[qname]][j,kk] <<- oldRow[kk]
    model$calculate(deps)
    return()
  }
  
  p <- w / s
  mstar <- as.integer(rcat(1, prob = p))

  for (kk in 1:K) model[[qname]][j,kk] <<- candidate[mstar, kk]
  model$calculate(deps)
},
methods = list(reset = function(){})
)

### 2. One simulation -------------------------------------------------------

run_one_simulation_nimble_unknownQ <- function(seed_number = 123,
                                               N = 200,
                                               J = 6,
                                               n.chains = 2,
                                               n.burnin = 1000,
                                               n.iter   = 2000,
                                               thin     = 1,
                                               debug = FALSE) {
  
  set.seed(seed_number)
  K <- 3
  P <- 6
  N <- 200
  J <- 18
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
      1,0,1,
      1,0,0,
      0,1,0,
      0,0,1
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
      
      0,1,1,
      1,0,0,
      1,0,0,
      1,0,0,
      0,1,0,
      0,0,1
    ), nrow = 18, byrow = TRUE)
  
  Q1_true <- Q1
  Q2_true <- Q2
  Q3_true <- Q3
  
 # theta_alpha <- 8
 # theta_beta  <- 12
   theta_alpha <- 1
   theta_beta  <- 29
  
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
   gamma01_12 <- matrix(
    c(
      -0.021, -0.681, -0.192,  0.083,  0.379, -0.104, -0.075,
      0.024, -0.054,  0.109,  0.191,  0.059,  0.125, -0.456,
      0.10,  -0.20,   0.05,   0.15,  -0.10,   0.05,  -0.30
    ),
    nrow = K, byrow = TRUE
  )
  gamma01_23 <- matrix(
    c(
      -0.121, -0.081, -0.092,  0.183,  0.379, -0.204, -0.175,
      0.224, -0.254,  0.069,  0.113,  -0.059,  0.125, -0.256,
      0.50,  -0.10,   -0.05,   0.35,  0.10,   0.05,  -0.10
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
  constants <- list(
    N = N, J = J, K = K, P = P
  )
  
  data_list <- list(
    Z = Z, Y1 = Y1, Y2 = Y2, Y3 = Y3
  )
  
  # params monitor 
  params_monitor <- c(
    "theta",
    "Q1_est","Q2_est","Q3_est",
    "g1","s1","g2","s2","g3","s3",
    "beta0","betaZ","gamma01_12","gamma01_23",
    "alpha1","alpha2","alpha3"
  )
  
  # Q needs to have two identity submatrices for initial values to ensure convergence
  Q_patterns <- rbind(
    c(1,0,0), # 1
    c(0,1,0), # 2
    c(0,0,1), # 3
    c(1,1,0), # 4
    c(1,0,1), # 5
    c(0,1,1), # 6
    c(1,1,1)  # 7
  )
  random_Q <- function(J){
    Q_patterns[sample(1:7, J, replace=TRUE), ]
  }

  # Initial values
  make_inits <- function() {
    
    theta_start <- runif(1, 0.2, 0.4)
    
    ## Q1
    Q1_init <- random_Q(J)
    
    Q1_init[1:3,] <- diag(3)
    Q1_init[13:15,] <- diag(3)
    
    ## Q2
    Q2_init <- random_Q(J)
    
    Q2_init[1:3,] <- diag(3)
    Q2_init[16:18,] <- diag(3)
    
    ## Q3
    Q3_init <- random_Q(J)
    
    Q3_init[1:3,] <- diag(3)
    Q3_init[7:9,] <- diag(3)
    
    # once mastery not forget thus alpha3 >= alpha2 >= alpha1
   # alpha1_init <- matrix(rbinom(N*K,1,0.5), N, K)
    
   # alpha2_init <- alpha1_init
    #flip <- matrix(rbinom(N*K,1,0.3), N, K)
   # alpha2_init[flip==1] <- 1
    
   # alpha3_init <- alpha2_init
  #  flip2 <- matrix(rbinom(N*K,1,0.3), N, K)
  #  alpha3_init[flip2==1] <- 1
    
    list(
      theta = theta_start,
      
      Q1_est = Q1_init,
      Q2_est = Q2_init,
      Q3_est = Q3_init,
      
    #  alpha1 = alpha1_init,
     # alpha2 = alpha2_init,
    #  alpha3 = alpha3_init,
      
      g1 = runif(J, 0.05, 0.1), s1 = runif(J, 0.05, 0.1),
      g2 = runif(J, 0.05, 0.1), s2 = runif(J, 0.05, 0.1),
      g3 = runif(J, 0.05, 0.1), s3 = runif(J, 0.05, 0.1),
      
      beta0 = rnorm(K, -1.1, 0.05),
      betaZ = matrix(rnorm(K * P, -0.2, 0.2), K, P),
      gamma01_12 = matrix(rnorm(K * (P + 1), 0, 0.1), K, (P + 1)),
      gamma01_23 = matrix(rnorm(K * (P + 1), 0, 0.1), K, (P + 1))
      
    )
  }
  
  # multiple chains initial values
  inits_list <- replicate(n.chains, make_inits(), simplify = FALSE)
  
  # build nimble model
  Rmodel <- nimbleModel(code = dina_lta_model,
                        constants = constants,
                        data = data_list,
                        inits = inits_list[[1]])  

  #Rmodel$calculate()
  
  # MCMC
  conf <- configureMCMC(Rmodel, monitors = params_monitor, thin = thin)
  # replace default samplers 
  for (j in 1:J) {
    for (k in 1:K) {
    conf$removeSampler(paste0("Q1_est[",j,",",k,"]"))
    conf$removeSampler(paste0("Q2_est[", j, ",", k, "]"))
    conf$removeSampler(paste0("Q3_est[", j, ",", k, "]"))
    }
  }
  
  # add row-wise sampler for each row
  for (j in 1:J) {
    conf$addSampler(
      target  = paste0("Q1_est[",j,",1:",K,"]"),
      type    = sampler_Q_rowwise,
      control = list(j=j, J=J, K=K, qname="Q1_est")
    )
    
    conf$addSampler(
      target  = paste0("Q2_est[",j,",1:",K,"]"),
      type    = sampler_Q_rowwise,
      control = list(j=j, J=J, K=K, qname="Q2_est")
    )
    
    conf$addSampler(
      target  = paste0("Q3_est[",j,",1:",K,"]"),
      type    = sampler_Q_rowwise,
      control = list(j=j, J=J, K=K, qname="Q3_est")
    )
  }
  
  # compile model: nimble package transfers nimbleCode and numbleFunction sampler to C++
  Cmodel <- compileNimble(Rmodel)
  
  # build MCMC
  Rmcmc <- buildMCMC(conf)
  # compile mcmc
  Cmcmc <- compileNimble(Rmcmc, project = Cmodel, showCompilerOutput = TRUE)
  
  debug <- isTRUE(debug)
  if (debug) {
    samp <- Cmcmc$mvSamples
    print(dim(as.matrix(samp)))
  }
  samp <- Cmcmc$mvSamples

  
  # run multiple MCMC chains
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
    success = TRUE,
    seed = seed_number,
    samples = samples,
    smat = smat,
    true = list(
      Q1_true = Q1_true, Q2_true = Q2_true, Q3_true = Q3_true,
      gs_true_t1 = gs_true_t1, gs_true_t2 = gs_true_t2, gs_true_t3 = gs_true_t3,
      beta0 = beta0, betaZ = betaZ, gamma01_12 = gamma01_12,  gamma01_23 = gamma01_23, 
      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3,
      N = N, J = J, K = K
    )
  )
}


### 3. Run model -------------------------------

out <- run_one_simulation_nimble_unknownQ(
  seed_number = 1,
  N = 200,
  J = 18,
  n.chains = 4,
  n.burnin = 1000,
  n.iter   = 1000,
  thin     = 5
)

str(out, max.level = 1)
J <- 18
K <- 3

### 4. Convergence diagnostics -------------------------------

diag_mcmc <- function(samples, rhat_thr = 1.10) {
  # 1) Rhat
  gd <- gelman.diag(samples, multivariate = FALSE)
  psrf_mat <- gd$psrf
  rhat <- psrf_mat[, 1]
  names(rhat) <- rownames(psrf_mat)
  
  # 2) ESS
  ess <- effectiveSize(samples)
  
  # 3) summary
  sum_df <- data.frame(
    max_rhat     = max(rhat, na.rm = TRUE),
    median_rhat  = median(rhat, na.rm = TRUE),
    num_rhat_bad = sum(rhat > rhat_thr, na.rm = TRUE),
    min_ess      = min(ess, na.rm = TRUE),
    median_ess   = median(ess, na.rm = TRUE),
    mean_ess     = mean(ess, na.rm = TRUE)
  )
  
  bad <- data.frame(
    parameter = names(rhat)[which(rhat > rhat_thr)],
    rhat = as.numeric(rhat[which(rhat > rhat_thr)])
  )
  
  list(summary = sum_df, bad_rhat = bad)
}

library(coda) # to use gelman.diag
diag_results <- diag_mcmc(out$samples) # this will take few seconds
diag_results$summary
diag_results$bad_rhat

# tract plots
par_trace <- c(
  "g1[1]","s1[1]","g2[1]","s2[1]","g3[1]","s3[1]")
traceplot(out$samples[, par_trace], main = "Trace plots: structural parameters")

q_trace <- c(
  "Q1_est[11, 1]","Q1_est[11, 2]","Q1_est[11, 3]",
  "Q1_est[12, 1]","Q1_est[12, 2]","Q1_est[12, 3]",
  "Q1_est[13, 1]","Q1_est[13, 2]","Q1_est[13, 3]"
)
traceplot(out$samples[, q_trace], main="Trace plots: selected Q1 entries")

### 5. Evaluation -------------------------------

## ---- 5a. help function -------------------------------------------------------
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

## PIP for Q (posterior inclusion probabilities)
get_Q_pip <- function(smat, varname, J, K){
  if(varname %in% c("Q1","Q2","Q3")) varname <- paste0(varname,"_est")
  cn <- colnames(smat)
  idx <- grep(paste0("^",varname,"\\["), cn)
  
  ij <- do.call(rbind,
                regmatches(cn[idx],
                           regexec("\\[(\\d+)\\s*,\\s*(\\d+)\\]", cn[idx])))
  
  pip <- matrix(NA_real_, J, K)
  pip[cbind(as.integer(ij[,2]), as.integer(ij[,3]))] <-
    colMeans(smat[,idx,drop=FALSE])
  pip
}

## Evaluation for Q (continuous + PIP)
eval_Q <- function(pip, Q_true, ambig_band = 0.05) {
  # pip: J×K posterior inclusion probabilities, [0,1]
  # continuous error on probabilities 
  diff <- as.numeric(pip) - as.numeric(Q_true)
  rmse <- sqrt(mean(diff^2))
  mae  <- mean(abs(diff))
  
  list(
    rmse = rmse,
    mae  = mae,
    mean_pip_on_ones  = mean(pip[Q_true == 1]),
    mean_pip_on_zeros = mean(pip[Q_true == 0]),
    prop_pip_ambig    = mean(abs(pip - 0.5) < ambig_band)
  )
}

## Binary Q_hat from PIP threshold
pip_to_Qhat <- function(pip, thr = 0.5) 1L * (pip >= thr)


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


compute_Q_metrics <- function(Q_hat, Q_true) {
  est  <- as.integer(Q_hat)
  true <- as.integer(Q_true)
  
  TP <- sum(est == 1 & true == 1)
  FP <- sum(est == 1 & true == 0)
  TN <- sum(est == 0 & true == 0)
  FN <- sum(est == 0 & true == 1)
  
  # avoid x/0
  safe_div <- function(a, b) if (b == 0) NA_real_ else a / b
  
  FPR <- safe_div(FP, FP + TN)
  FNR <- safe_div(FN, TP + FN)
  ACC <- safe_div(TP + TN, TP + TN + FP + FN)
  
  row_ACC <- mean(apply(Q_hat == Q_true, 1, all)) 
  
  data.frame(
    FPR = FPR,
    FNR = FNR,
    ACC = ACC,
    row_ACC = row_ACC
  )
}

## Prior for theta
## candidates in fixed order:
## 1:001 2:010 3:011 4:100 5:101 6:110 7:111
cand_rowSums <- function(K = 3) c(1,1,2,1,2,2,3)

# Marginal prior: integrate over theta ~ Beta(a,b)
# p(q) ~ B(a+s, b+K-s) / B(a,b)
prior_probs_marginal <- function(a, b, K = 3){
  s <- cand_rowSums(K)
  w <- beta(a + s, b + K - s) / beta(a, b)
  w / sum(w)
}

# Conditional prior: given theta
# p(q|theta) ~ theta^s (1-theta)^(K-s)
# --------------------------
prior_probs_conditional <- function(theta, K = 3){
  s <- cand_rowSums(K)
  w <- theta^s * (1 - theta)^(K - s)
  w / sum(w)
}

# examples
#prior_probs_marginal(a=8, b=12, K=3)
#prior_probs_marginal(a=1, b=29, K=3)

#prior_probs_conditional(theta=0.4, K=3)
#prior_probs_conditional(theta=0.1, K=3)


## ---- 5b. main evaluation function --------------------------------------
evaluate_one_run <- function(smat,
                             J, K, N,
                             Q1_true, Q2_true, Q3_true,
                             gs_true_t1, gs_true_t2, gs_true_t3,
                             beta0, betaZ, gamma01_12,gamma01_23,
                             alpha1, alpha2, alpha3) {
  
  ## 1) Q: PIP + binary Q_hat
  Q1_pip <- get_Q_pip(smat, "Q1", J, K)
  Q2_pip <- get_Q_pip(smat, "Q2", J, K)
  Q3_pip <- get_Q_pip(smat, "Q3", J, K)
  
  Q1_hat <- pip_to_Qhat(Q1_pip, thr = 0.5)
  Q2_hat <- pip_to_Qhat(Q2_pip, thr = 0.5)
  Q3_hat <- pip_to_Qhat(Q3_pip, thr = 0.5)
  
  ## 2) Item parameters (posterior means)
  g1_est <- get_posterior_mean_vec(smat, "^g1\\[")
  s1_est <- get_posterior_mean_vec(smat, "^s1\\[")
  g2_est <- get_posterior_mean_vec(smat, "^g2\\[")
  s2_est <- get_posterior_mean_vec(smat, "^s2\\[")
  g3_est <- get_posterior_mean_vec(smat, "^g3\\[")
  s3_est <- get_posterior_mean_vec(smat, "^s3\\[")
  
  ## 3) Regression parameters
  beta0_est <- get_posterior_mean_vec(smat, "^beta0\\[")          # length K
  betaZ_est <- get_posterior_mean_mat(smat, "^betaZ\\[", K, 6)    # K x P (P=6)
  gamma01_12_est <- get_posterior_mean_mat(smat, "^gamma01_12\\[", K, 7)# K x (P+1)
  gamma01_23_est <- get_posterior_mean_mat(smat, "^gamma01_23\\[", K, 7)# K x (P+1)
  
  ## 4) Latent attributes: posterior prob + hard classification
  a1_postprob <- get_posterior_mean_mat(smat, "^alpha1\\[", N, K)
  a2_postprob <- get_posterior_mean_mat(smat, "^alpha2\\[", N, K)
  a3_postprob <- get_posterior_mean_mat(smat, "^alpha3\\[", N, K)
  a1_hat <- 1L * (a1_postprob >= 0.5)
  a2_hat <- 1L * (a2_postprob >= 0.5)
  a3_hat <- 1L * (a3_postprob >= 0.5)
  
  ## 5) Recovery metrics
  ## 5a) Q recovery
  # binary + PIP diagnostics (summarize the ones and zeros of each entry of Q in the posterior samples after discarding the burn-in samples)
  res_Q1 <- eval_Q(Q1_pip, Q1_true)
  res_Q2 <- eval_Q(Q2_pip, Q2_true)
  res_Q3 <- eval_Q(Q3_pip, Q3_true)
  
  # Confusion metrics (TPR/FPR/FNR/ACC)
  ACC_Q1 <- compute_Q_metrics(Q1_hat, Q1_true)
  ACC_Q2 <- compute_Q_metrics(Q2_hat, Q2_true)
  ACC_Q3 <- compute_Q_metrics(Q3_hat, Q3_true)
  
  # Q_results_recovery list
  Q1_res_recovery <- list(
    continuous_pip = res_Q1, metrics = ACC_Q1)
  Q2_res_recovery <- list(
    continuous_pip = res_Q2, metrics = ACC_Q2)
  Q3_res_recovery <- list(
    continuous_pip = res_Q3, metrics = ACC_Q3)
  
  ## 5b) item parameters
  res_g1 <- eval_cont_vec(g1_est, gs_true_t1[, 1])
  res_s1 <- eval_cont_vec(s1_est, gs_true_t1[, 2])
  res_g2 <- eval_cont_vec(g2_est, gs_true_t2[, 1])
  res_s2 <- eval_cont_vec(s2_est, gs_true_t2[, 2])
  res_g3 <- eval_cont_vec(g3_est, gs_true_t3[, 1])
  res_s3 <- eval_cont_vec(s3_est, gs_true_t3[, 2])
  
  ## 5c) regression parameters
  res_beta0   <- eval_cont_vec(beta0_est, beta0)
  res_betaZ   <- eval_cont_mat(betaZ_est, betaZ)
  res_gamma01_12 <- eval_cont_mat(gamma01_12_est, gamma01_12)
  res_gamma01_23 <- eval_cont_mat(gamma01_23_est, gamma01_23)
  
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
  # a1_postprob[17,]
  
  ## 6) Uncertainty: 95% CIs for key continuous parameters
  ci <- list(
    theta = get_ci_mat(smat, "^theta"),
    g1 = get_ci_mat(smat, "^g1\\["),
    s1 = get_ci_mat(smat, "^s1\\["),
    g2 = get_ci_mat(smat, "^g2\\["),
    s2 = get_ci_mat(smat, "^s2\\["),
    g3 = get_ci_mat(smat, "^g3\\["),
    s3 = get_ci_mat(smat, "^s3\\["),
    beta0  = get_ci_mat(smat, "^beta0\\["),
    betaZ  = get_ci_mat(smat, "^betaZ\\["),
    gamma01_12 = get_ci_mat(smat, "^gamma01_12\\["),
    gamma01_23 = get_ci_mat(smat, "^gamma01_23\\[")
  )
  
  ## 7) prior for theta
  theta_mean <- mean(smat[,"theta"])
  
  prior_theta <- list(
    marginal = prior_probs_marginal(a=1, b=29, K=K),
    conditional_at_theta_mean = prior_probs_conditional(theta_mean, K=K),
    theta_mean = theta_mean
  )
  
  list(
    
    ## 1) estimates for key parameters
    est = list(
      ## Q
      Q1_pip = Q1_pip, Q2_pip = Q2_pip, Q3_pip = Q3_pip,
      Q1_hat = Q1_hat, Q2_hat = Q2_hat, Q3_hat = Q3_hat,
      
      ## item parameters
      g1 = g1_est, s1 = s1_est,
      g2 = g2_est, s2 = s2_est,
      g3 = g3_est, s3 = s3_est,
      
      ## regressions
      beta0 = beta0_est,
      betaZ = betaZ_est,
      gamma01_12 = gamma01_12_est,
      gamma01_23 = gamma01_23_est#,
      
      ## latent attributes
   #   a1_postprob = a1_postprob, a2_postprob = a2_postprob, a3_postprob = a3_postprob,
   #   a1_hat = a1_hat, a2_hat = a2_hat, a3_hat = a3_hat
    ),
    
    ## 2) recovery / evaluation metrics
    metrics = list(
      ## Q: binary + pip diagnostics
     Q_binpip = list(Q1 = res_Q1, Q2 = res_Q2, Q3 = res_Q3),
      
      ## Q: confusion-style metrics (TPR/FPR/FNR/ACC)
      Q_confusion = list(Q1 = ACC_Q1, Q2 = ACC_Q2, Q3 = ACC_Q3),
      
      ## Q: bundled (if you like this format)
      Q_recovery = list(
        Q1 = Q1_res_recovery,
        Q2 = Q2_res_recovery,
        Q3 = Q3_res_recovery
      ),
      
      ## item parameters
      g1 = res_g1, s1 = res_s1,
      g2 = res_g2, s2 = res_s2,
      g3 = res_g3, s3 = res_s3,
      
      ## regression parameters
      beta0 = res_beta0,
      betaZ = res_betaZ,
      gamma01_12 = res_gamma01_12,
      gamma01_23 = res_gamma01_23,
      
      ## attribute classification
      PAR = c(t1 = PAR_t1, t2 = PAR_t2, t3 = PAR_t3),
      AAR = rbind(t1 = AAR_t1, t2 = AAR_t2, t3 = AAR_t3)#,
      
      ## individual attribute profile diagnostics (select few to check)
  #    alpha_diag = list(t1 = diag1, t2 = diag2, t3 = diag3)
    ),
    
    ## 3) uncertainty summaries (95% CIs)
    ci = ci,
    ## 4) prior for theta 
    prior_theta
  )
}


res1 <- evaluate_one_run(
  smat = out$smat,
  J = out$true$J, K = out$true$K, N = out$true$N,
  Q1_true = out$true$Q1_true, Q2_true = out$true$Q2_true, Q3_true = out$true$Q3_true,
  gs_true_t1 = out$true$gs_true_t1, gs_true_t2 = out$true$gs_true_t2, gs_true_t3 = out$true$gs_true_t3,
  beta0 = out$true$beta0, betaZ = out$true$betaZ, gamma01_12 = out$true$gamma01_12, gamma01_23 = out$true$gamma01_23, 
  alpha1 = out$true$alpha1, alpha2 = out$true$alpha2, alpha3 = out$true$alpha3
)

res1


### 8. Convergence plots---------------------
library(coda)
running_mean_plot <- function(samples, param, main=NULL){
  # samples: mcmc.list
  mlist <- samples
  par(mfrow=c(length(mlist),1), mar=c(3,4,2,1))
  for(ch in seq_along(mlist)){
    x <- as.numeric(mlist[[ch]][, param])
    rm <- cumsum(x) / seq_along(x)
    plot(rm, type="l",
         xlab="Iteration (post-thin index)",
         ylab="Running mean",
         main = if(is.null(main)) paste0("Chain ", ch, ": ", param) else main)
    abline(h = mean(x), lty=2)
  }
  par(mfrow=c(1,1))
}

running_mean_plot(out$samples, "theta")
running_mean_plot(out$samples, "beta0[1]")
running_mean_plot(out$samples, "betaZ[1, 1]")


# Gelman-Rubin iteration plots
library(coda)

pars_core <- c(
  "theta",
  paste0("beta0[",1:3,"]"),
  "g1[1]","s1[1]","g2[1]","s2[1]","g3[1]","s3[1]"
)

gelman.plot(out$samples[, pars_core], autoburnin = FALSE)
gelman.diag(out$samples[, pars_core], multivariate = FALSE)
effectiveSize(out$samples[, pars_core])




