# Simulation Example

### 0. Load Packages -----------------------------------------------------------
library(rjags)   # For MCMC sampling using JAGS
library(coda)    # For convergence diagnostics
library(MASS)    # For multivariate normal distribution (mvrnorm)
library(dplyr)   # For data cleaning
library(pROC)    # For computing AUC
library(tidyr)   # For data reshaping

jags_model_string_unknownQ <- "
model {

  #### Hyperpriors for Q-matrix sparsity
  theta ~ dbeta(theta_alpha, theta_beta) T(1e-6, 1 - 1e-6)

  #### Q-matrix estimation (free from j=3 to J; first two rows fixed for identifiability)
  for (j in 3:J) {
    for (k in 1:K) {
      Q1_est[j,k] ~ dbeta(theta + 0.001, (1 - theta) + 0.001) T(1e-6, 1 - 1e-6)
      Q2_est[j,k] ~ dbeta(theta + 0.001, (1 - theta) + 0.001) T(1e-6, 1 - 1e-6)
    }
  }

  # Fixed rows in Q-matrix for identifiability
  Q1_est[1,1] <- 1; Q1_est[1,2] <- 0
  Q1_est[2,1] <- 0; Q1_est[2,2] <- 1
  Q2_est[1,1] <- 0; Q2_est[1,2] <- 1
  Q2_est[2,1] <- 1; Q2_est[2,2] <- 0

  #### Item parameters (guessing and slipping)
  for (j in 1:J) {
    g1[j] ~ dbeta(1,1)
    s1[j] ~ dbeta(1,1)
    g2[j] ~ dbeta(1,1)
    s2[j] ~ dbeta(1,1)
  }

  #### Regression coefficients for latent attributes
  for (k in 1:K) {
    beta0[k] ~ dnorm(0, 1)
    for (p in 1:P) {
      betaZ[k,p] ~ dnorm(0, 1)
    }
    gamma01[k,1] ~ dnorm(0, 1)
    gamma10[k,1] ~ dnorm(0, 1)
    for (p in 2:(P + 1)) {
      gamma01[k,p] ~ dnorm(0, 1)
      gamma10[k,p] ~ dnorm(0, 1)
    }
  }

  #### Latent attributes at time 1 and transition to time 2
  for (n in 1:N) {
    for (k in 1:K) {
      logit(p1[n,k]) <- beta0[k] + inprod(betaZ[k,], Z[n,])
      alpha1[n,k] ~ dbern(p1[n,k])
    }

    for (k in 1:K) {
      logit(p01[n,k]) <- gamma01[k,1] + inprod(gamma01[k,2:(P+1)], Z[n,])
      logit(p10[n,k]) <- gamma10[k,1] + inprod(gamma10[k,2:(P+1)], Z[n,])
      prob_alpha2[n,k] <- alpha1[n,k] * (1 - p10[n,k]) + (1 - alpha1[n,k]) * p01[n,k]
      alpha2[n,k] ~ dbern(prob_alpha2[n,k])
    }
  }

  #### Measurement model for item responses at both time points
  for (n in 1:N) {
    for (j in 1:J) {
      for (k in 1:K) {
        tmp_log1[n,j,k] <- Q1_est[j,k] * log(alpha1[n,k] + 1e-12)
      }
      log_eta1[n,j] <- sum(tmp_log1[n,j,])
      eta1[n,j] <- exp(log_eta1[n,j])
      pY1[n,j] <- (1 - s1[j]) * eta1[n,j] + g1[j] * (1 - eta1[n,j])
      Y1[n,j] ~ dbern(pY1[n,j])

      for (k in 1:K) {
        tmp_log2[n,j,k] <- Q2_est[j,k] * log(alpha2[n,k] + 1e-12)
      }
      log_eta2[n,j] <- sum(tmp_log2[n,j,])
      eta2[n,j] <- exp(log_eta2[n,j])
      pY2[n,j] <- (1 - s2[j]) * eta2[n,j] + g2[j] * (1 - eta2[n,j])
      Y2[n,j] ~ dbern(pY2[n,j])
    }
  }

  #### Save posterior probabilities
  for (n in 1:N) {
    for (k in 1:K) {
      prob_attr1[n,k] <- p1[n,k]
      prob_attr2[n,k] <- prob_alpha2[n,k]
    }
  }
}"



### 2. Helpful functions ---------------------------------------------------------
# Evaluate RMSE and bias between estimated and true continuous matrices
evaluate_continuous_matrix <- function(est, true) {
  rmse <- sqrt(mean((est - true)^2))
  bias <- mean(est - true)
  list(rmse = rmse, bias = bias)
}

# Evaluate binary matrix estimation accuracy and RMSE
evaluate_binary_matrix <- function(est, true) {
  accuracy <- mean(round(est) == true)
  rmse     <- sqrt(mean((est - true)^2))
  list(accuracy = accuracy, rmse = rmse)
}

# Compute confusion matrix metrics: TP, FP, TN, FN, FPR, FNR, TPR, TNR, ACC
compute_confusion <- function(est, true) {
  v_est  <- as.integer(est > 0.5)
  v_true <- as.integer(true)
  TP <- sum(v_est &  v_true)
  FP <- sum(v_est & !v_true)
  TN <- sum(!v_est & !v_true)
  FN <- sum(!v_est &  v_true)
  data.frame(TP, FP, TN, FN,
             FPR = FP / (FP + TN),
             FNR = FN / (FN + TP),
             TPR = TP / (TP + FN),
             TNR = TN / (TN + FP),
             ACC = (TP + TN) / (TP + FP + TN + FN))
}
# Compute AUC for each skill dimension
get_auc_vec <- function(prob_mat, truth_mat, prefix = "t") {
  K  <- ncol(truth_mat)
  au <- numeric(K)
  for (k in seq_len(K)) {
    y <- truth_mat[, k]
    p <- prob_mat[, k]
    keep <- !(is.na(y) | is.na(p))
    y <- y[keep]; p <- p[keep]
    if (length(y) == 0L || length(unique(y)) < 2L || all(p == p[1L])) {
      au[k] <- NA
      next
    }
    au[k] <- tryCatch({
      as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
    }, error = function(e) NA_real_)
  }
  names(au) <- sprintf("%s_k%d", prefix, seq_len(K))
  au
}
# Generate fixed Q-matrix pairs for time 1 and time 2 based on total items J
make_Q_pair_fixed <- function(J) {
  stopifnot(J %in% c(6, 18, 30))
  K <- 2
  if (J == 6) {
    Q1 <- rbind(diag(1, K), diag(1, K), c(1, 0), c(1, 1))
    Q2 <- rbind(c(0, 1), c(1, 0), c(0, 1), c(1, 0), c(0, 1), c(1, 1))
  }
  if (J == 18) {
    Q1 <- rbind(matrix(rep(diag(1, K), 3), byrow = TRUE, ncol = K),
                matrix(rep(diag(1, K), 3), byrow = TRUE, ncol = K),
                c(1, 0), c(0, 1),
                matrix(rep(c(1, 1), 4), ncol = K, byrow = TRUE))
    Q2 <- rbind(c(0, 1), c(1, 0),
                matrix(rep(c(0, 1), 6), ncol = K, byrow = TRUE),
                matrix(rep(c(1, 0), 6), ncol = K, byrow = TRUE),
                matrix(rep(c(1, 1), 4), ncol = K, byrow = TRUE))
  }
  if (J == 30) {
    Q1 <- rbind(matrix(rep(diag(1, K), 6), ncol = K, byrow = TRUE),
                matrix(rep(diag(1, K), 6), ncol = K, byrow = TRUE),
                matrix(rep(c(1, 1), 6), ncol = K, byrow = TRUE))
    Q2 <- rbind(c(0, 1), c(1, 0),
                matrix(rep(c(0, 1), 11), ncol = K, byrow = TRUE),
                matrix(rep(c(1, 0), 11), ncol = K, byrow = TRUE),
                matrix(rep(c(1, 1), 6),  ncol = K, byrow = TRUE))
  }
  list(Q_time1 = Q1, Q_time2 = Q2)
}



### 
run_one_simulation_rjags_unknownQ <- function(seed_number = 123,
                                              N = 200, 
                                              J = 6,
                                              n.chains = 4,
                                              n.adapt  = 1000,
                                              n.burnin = 1000,
                                              n.iter   = 2000,
                                              thin     = 1) {
  set.seed(seed_number)
  K <- 2;  P <- 6
  
  ## ----------------------------------------------
  ## (i) Data Generation
  ## ----------------------------------------------
  Z <- mvrnorm(N, mu = rep(0, P), Sigma = diag(P))
  Q_pair <- make_Q_pair_fixed(J)
  Q_t1_true <- Q_pair$Q_time1
  Q_t2_true <- Q_pair$Q_time2
  
  Q_true <- Q_pair$Q_time1   # or Q_time2
  rho <- sum(Q_true != 0) / ((J - 2) * K)  # minus the fixed two rows
  strength <- 20
  theta_alpha <- rho * strength
  theta_beta  <- (1 - rho) * strength
  
  # True guessing and slipping parameters (for both time points)
  gs_true_t1 <- matrix(runif(J * 2, 0.05, 0.20), ncol = 2)
  gs_true_t2 <- matrix(runif(J * 2, 0.05, 0.20), ncol = 2)
  
  # True regression coefficients (beta0 and betaZ)
  beta0 <- c(-1.12, -1.193)
  betaZ <- matrix(c(
    -0.576,  0.092,  1.282, -0.03,  -0.196,  0.06,
    0.125, -0.302, -0.03,   1.315, -0.567,  0.261
  ), nrow = K, byrow = TRUE)
  
  # Transition parameters: gamma01 (0->1), gamma10 (1->0)
  gamma01 <- matrix(c(
    -0.021, -0.681, -0.192,  0.083,  0.379, -0.104, -0.075,
    0.024, -0.054,  0.109,  0.191,  0.059,  0.125, -0.456
  ), nrow = K, byrow = TRUE)
  
  gamma10 <- matrix(c(
    -0.186, -0.241, -0.511, -0.10, -0.189,  0.210, -0.384,
    -0.142, -0.165, -0.311, -0.071, -0.471, -0.155, -0.208
  ), nrow = K, byrow = TRUE)
  
  # Simulate latent attributes at t1 and t2
  invlogit <- function(x) 1 / (1 + exp(-x))
  p_init <- t(apply(Z, 1, function(z) invlogit(beta0 + betaZ %*% z)))
  alpha1 <- matrix(rbinom(N * K, 1, as.vector(p_init)), N, K)
  alpha2 <- matrix(0, N, K)
  for (i in 1:N) {
    z_ext <- c(1, Z[i, ])
    for (k in 1:K) {
      p01 <- invlogit(sum(gamma01[k, ] * z_ext))
      p10 <- invlogit(sum(gamma10[k, ] * z_ext))
      alpha2[i, k] <- if (alpha1[i, k] == 0) rbinom(1, 1, p01) else rbinom(1, 1, 1 - p10)
    }
  }
  
  # Generate item responses at t1 and t2
  gen_Y <- function(alpha, Q, gs) {
    outer(seq_len(N), seq_len(nrow(Q)), Vectorize(function(i, j) {
      eta <- prod(alpha[i, ] ^ Q[j, ])
      p   <- (1 - gs[j, 2]) * eta + gs[j, 1] * (1 - eta)
      rbinom(1, 1, p)
    }))
  }
  Y1 <- gen_Y(alpha1, Q_pair$Q_time1, gs_true_t1)
  Y2 <- gen_Y(alpha2, Q_pair$Q_time2, gs_true_t2)
  
  ## ----------------------------------------------
  ## (ii) Fit JAGS Model
  ## ----------------------------------------------
  
  jags_data <- list(
    N = N, J = J, K = K, P = P,
    Z = Z, Y1 = Y1, Y2 = Y2,
    theta_alpha = theta_alpha,
    theta_beta = theta_beta
  )
  
  params_monitor <- c("theta","Q1_est","Q2_est",
                      "g1", "s1", "g2", "s2",
                      "beta0", "betaZ", "gamma01", "gamma10",
                      "alpha1", "alpha2")
  # Function to generate initial values
  make_inits <- function() { 
    theta_start <- runif(1, 0.2, 0.8)  
    epsilon <- 0.05  
    Q1_init <- matrix(rbeta(J*K, theta_start + epsilon, (1 - theta_start) + epsilon), J, K)
    Q2_init <- matrix(rbeta(J*K, theta_start + epsilon, (1 - theta_start) + epsilon), J, K)
    Q1_init[1:2,] <- NA  
    Q2_init[1:2,] <- NA 
    
    list(
      theta   = theta_start,
      Q1_est  = Q1_init,
      Q2_est  = Q2_init,
      g1      = runif(J,0,0.3), s1 = runif(J,0,0.3),
      g2      = runif(J,0,0.3), s2 = runif(J,0,0.3),
      beta0   = c(-1.1, -1.2) + rnorm(K, 0, 0.1),
      gamma01 = gamma01 + matrix(rnorm(K*(P+1),0,0.1), K, P+1),
      gamma10 = gamma10 + matrix(rnorm(K*(P+1),0,0.1), K, P+1)
    )
  }
  
  # Run JAGS model
  jm <- jags.model(
    textConnection(jags_model_string_unknownQ),
    data     = jags_data,
    inits    = replicate(n.chains, make_inits(), FALSE),
    n.chains = n.chains,
    n.adapt  = n.adapt
  )
  
  update(jm, n.burnin)
  
  samples <- coda.samples(
    jm, params_monitor,
    n.iter = n.iter, thin = thin
  )
  
  
  ## ----------------------------------------------
  ## (iii) Extract Parameter Estimates
  ## ----------------------------------------------
  smat    <- as.matrix(samples)
  col_mean <- colMeans(smat)
  
  ## ---- Posterior means for Q
  Q1_cols <- grep("^Q1_est\\[", names(col_mean))
  Q2_cols <- grep("^Q2_est\\[", names(col_mean))
  Q_est_t1 <- matrix(col_mean[Q1_cols], nrow=J, ncol=K, byrow=FALSE)
  Q_est_t2 <- matrix(col_mean[Q2_cols], nrow=J, ncol=K, byrow=FALSE)
  
  ## ---- Posterior means for items
  g1_est <- col_mean[grep("^g1\\[", names(col_mean))]
  s1_est <- col_mean[grep("^s1\\[", names(col_mean))]
  g2_est <- col_mean[grep("^g2\\[", names(col_mean))]
  s2_est <- col_mean[grep("^s2\\[", names(col_mean))]
  
  ## ---- Posterior means for regression
  gamma01_est <- col_mean[grep("^gamma01\\[", names(col_mean))]
  gamma10_est <- col_mean[grep("^gamma10\\[", names(col_mean))]
  gamma01_est_mat <- matrix(gamma01_est, nrow=K, byrow=FALSE)
  gamma10_est_mat <- matrix(gamma10_est, nrow=K, byrow=FALSE)
  
  beta0_est <- col_mean[grep("^beta0\\[", names(col_mean))]
  betaZ_est <- col_mean[grep("^betaZ\\[", names(col_mean))]
  betaZ_est_mat <- matrix(betaZ_est, nrow=K, byrow=FALSE)
  
  ## ---- Posterior means for attributes
  pr1_cols <- grep("^alpha1\\[", names(col_mean))
  pr2_cols <- grep("^alpha2\\[", names(col_mean))
  prob_attr1 <- matrix(col_mean[pr1_cols], N, K)
  prob_attr2 <- matrix(col_mean[pr2_cols], N, K)
  alpha1_pred <- ifelse(prob_attr1 > 0.5, 1, 0)
  alpha2_pred <- ifelse(prob_attr2 > 0.5, 1, 0)
  
  # Posterior probability of alpha
  a1_cols <- grep("^alpha1\\[", colnames(smat))
  a2_cols <- grep("^alpha2\\[", colnames(smat))
  # Posterior mean (mode) for alpha
  alpha1_postprob <- matrix(colMeans(smat[, a1_cols]), N, K)
  alpha2_postprob <- matrix(colMeans(smat[, a2_cols]), N, K)
  
  ## ----------------------------------------------
  ## (iv) Evaluation Metrics
  ## ----------------------------------------------
  # Q-matrix recovery
  res_q1 <- evaluate_binary_matrix(Q_est_t1, Q_t1_true)
  res_q2 <- evaluate_binary_matrix(Q_est_t2, Q_t2_true)
  
  # Guessing and slipping parameter recovery
  res_g_t1 <- evaluate_continuous_matrix(col_mean[grep("^g1\\[", names(col_mean))], gs_true_t1[,1])
  res_s_t1 <- evaluate_continuous_matrix(col_mean[grep("^s1\\[", names(col_mean))], gs_true_t1[,2])
  res_g_t2 <- evaluate_continuous_matrix(col_mean[grep("^g2\\[", names(col_mean))], gs_true_t2[,1])
  res_s_t2 <- evaluate_continuous_matrix(col_mean[grep("^s2\\[", names(col_mean))], gs_true_t2[,2])
  
  # Transition parameters
  res_gamma01_K1 <- evaluate_continuous_matrix(gamma01_est_mat[1,], gamma01[1,])
  res_gamma01_K2 <- evaluate_continuous_matrix(gamma01_est_mat[2,], gamma01[2,])
  res_gamma10_K1 <- evaluate_continuous_matrix(gamma10_est_mat[1,], gamma10[1,])
  res_gamma10_K2 <- evaluate_continuous_matrix(gamma10_est_mat[2,], gamma10[2,])
  
  # Regression parameters
  res_beta0_K1 <- evaluate_continuous_matrix(beta0_est[1], beta0[1])
  res_beta0_K2 <- evaluate_continuous_matrix(beta0_est[2], beta0[2])
  res_betaZ_K1 <- evaluate_continuous_matrix(betaZ_est_mat[1,], betaZ[1,])
  res_betaZ_K2 <- evaluate_continuous_matrix(betaZ_est_mat[2,], betaZ[2,])

  # Attribute classification accuracy
  compute_PAR <- function(pred, true) mean(apply(pred == true, 1, all))
  compute_AAR <- function(pred, true) colMeans(pred == true)
  PAR_time1 <- compute_PAR(alpha1_pred, alpha1)
  PAR_time2 <- compute_PAR(alpha2_pred, alpha2)
  AAR_time1 <- compute_AAR(alpha1_pred, alpha1)
  AAR_time2 <- compute_AAR(alpha2_pred, alpha2)
  
  # AUC
  auc_t1 <- get_auc_vec(prob_attr1, alpha1, "t1")
  auc_t2 <- get_auc_vec(prob_attr2, alpha2, "t2")
  
  auc_t1_mode <- get_auc_vec(alpha1_postprob, alpha1, "t1")
  auc_t2_mode <- get_auc_vec(alpha2_postprob, alpha2, "t2")
  
  # Confusion metrics for alpha and Q
  cm_t1_k1 <- compute_confusion(alpha1_pred[,1], alpha1[,1])
  cm_t1_k2 <- compute_confusion(alpha1_pred[,2], alpha1[,2])
  cm_t2_k1 <- compute_confusion(alpha2_pred[,1], alpha2[,1])
  cm_t2_k2 <- compute_confusion(alpha2_pred[,2], alpha2[,2])
  
  flatten_cm <- function(df,prefix){
    metrics <- c("FPR","FNR","TPR","TNR","ACC")
    setNames(as.numeric(df[1,metrics]),
             paste0("alpha_metrics_",prefix,"_",metrics))
  }
  flat_alpha <- c(flatten_cm(cm_t1_k1,"t1_k1"),
                  flatten_cm(cm_t1_k2,"t1_k2"),
                  flatten_cm(cm_t2_k1,"t2_k1"),
                  flatten_cm(cm_t2_k2,"t2_k2"))
  
  # Confusion matrices for Q ----------
  mask_t1 <- !(1:J %in% c(1,2)) 
  mask_t2 <- !(1:J %in% c(1,2))
  
  cm_q_time1 <- compute_confusion(Q_est_t1[mask_t1,], Q_t1_true[mask_t1,])  # t = 1
  cm_q_time2 <- compute_confusion(Q_est_t2[mask_t2,], Q_t2_true[mask_t2,])  # t = 2
  
  # Organize the results from confusion matrices
  flatten_cm <- function(df, prefix){
    metrics <- c("FPR", "FNR", "TPR", "TNR", "ACC")
    setNames(as.numeric(df[1, metrics]),
             paste0("q_metrics_", prefix, "_", metrics))
  }
  
  flat_q <- c(flatten_cm(cm_q_time1, "t1"),
              flatten_cm(cm_q_time2, "t2"))
  
  
  ## ----------------------------------------------
  ## (v) Convergence Diagnostics 
  ## ----------------------------------------------
  # Compute R-hat (Gelman-Rubin statistic) for all parameters
  psrf     <- tryCatch(gelman.diag(samples, multivariate = FALSE)$psrf[, 1], error = function(e) NA)
  # Remove fixed Q-matrix rows (they may have R-hat = 0 or NA)
  is_fixedQ <- grepl("^Q[12]\\[", names(psrf)) & (psrf == 0 | is.nan(psrf))
  psrf  <- psrf[!is_fixedQ]
  # Identify parameters with poor convergence (R-hat > threshold)
  thr_rhat  <- 1.10
  idx_bad   <- which( psrf > thr_rhat | is.na(psrf))
  
  # Save problematic parameters
  if (length(idx_bad) > 0) {
    bad_params <- data.frame(
      parameter = names(psrf)[idx_bad],
      rhat      = psrf[idx_bad]
    )
    # print(bad_params)                
    # write.csv(bad_params, "rhat_over_1.10.csv", row.names = FALSE) # if needed
  } else {
    message("All chains converged (R-hat ≤ ", thr_rhat, ").")
  }
  
  max_rhat <- max(psrf, na.rm = TRUE)
  
  # Compute effective sample size (ESS)
  ess <- tryCatch(effectiveSize(samples), error = function(e) NA)
  ess <- ess[ess > 0] 
  
  # Summarize ESS statistics
  min_ess    <- min(ess, na.rm = TRUE)
  median_ess <- median(ess, na.rm = TRUE)
  mean_ess   <- mean(ess, na.rm = TRUE)
  max_ess    <- max(ess, na.rm = TRUE)
  
  
  psrf_all <- tryCatch(gelman.diag(samples, multivariate = FALSE)$psrf[, 1], error = function(e) NA)
  non_alpha_idx <- !grepl("^alpha[12]\\[", names(psrf_all))
  psrf_non_alpha <- psrf_all[non_alpha_idx]
  max_rhat_non_alpha    <- max(psrf_non_alpha, na.rm = TRUE)
  median_rhat_non_alpha <- median(psrf_non_alpha, na.rm = TRUE)
  mean_rhat_non_alpha   <- mean(psrf_non_alpha, na.rm = TRUE)
  
  # Extract parameters with R-hat > 1.10
  high_rhat_params <- names(psrf_non_alpha)[psrf_non_alpha > 1.10]
  high_rhat_vals   <- psrf_non_alpha[psrf_non_alpha > 1.10]
  rhat_exceed_df <- data.frame(parameter = high_rhat_params,
                               rhat      = high_rhat_vals)
  
  psrf_all <- tryCatch(
    gelman.diag(samples, multivariate = FALSE)$psrf[, 1],
    error = function(e) NA
  )
  
  rhat_df <- data.frame(
    parameter = names(psrf),
    rhat      = psrf_all
  )
  
  # Summary table for convergence diagnostics
  diag_df <- data.frame(
    seed_used  = seed_number,
    max_rhat_all         = max(psrf_all, na.rm = TRUE),
    max_rhat_non_alpha   = max_rhat_non_alpha,
    median_rhat_non_alpha = median_rhat_non_alpha,
    mean_rhat_non_alpha  = mean_rhat_non_alpha,
    num_rhat_above_1.10  = sum(psrf_non_alpha > 1.10, na.rm = TRUE),
    num_rhat_above_1.20  = sum(psrf_non_alpha > 1.20, na.rm = TRUE),
    min_ess    = min_ess,
    median_ess = median_ess,
    mean_ess   = mean_ess,
    max_ess    = max_ess
  )
  
  ## ----------------------------------------------
  ## (iv) Return Summary and Full Results
  ## ----------------------------------------------
  summary_df <- data.frame(
    seed_used = seed_number,
    N_used    = N,
    J_used    = J,
    Q_time1_acc  = res_q1$accuracy,
    Q_time1_rmse = res_q1$rmse,
    Q_time2_acc  = res_q2$accuracy,
    Q_time2_rmse = res_q2$rmse,
    
    g_t1_rmse = res_g_t1$rmse,  g_t1_bias = res_g_t1$bias,
    g_t2_rmse = res_g_t2$rmse,  g_t2_bias = res_g_t2$bias,
    s_t1_rmse = res_s_t1$rmse,  s_t1_bias = res_s_t1$bias,
    s_t2_rmse = res_s_t2$rmse,  s_t2_bias = res_s_t2$bias,
    
    gamma01_K1_rmse = res_gamma01_K1$rmse, gamma01_K1_bias = res_gamma01_K1$bias,
    gamma01_K2_rmse = res_gamma01_K2$rmse, gamma01_K2_bias = res_gamma01_K2$bias,
    gamma10_K1_rmse = res_gamma10_K1$rmse, gamma10_K1_bias = res_gamma10_K1$bias,
    gamma10_K2_rmse = res_gamma10_K2$rmse, gamma10_K2_bias = res_gamma10_K2$bias,
    
    beta0_K1_rmse = res_beta0_K1$rmse, beta0_K1_bias = res_beta0_K1$bias,
    beta0_K2_rmse = res_beta0_K2$rmse, beta0_K2_bias = res_beta0_K2$bias,
    betaZ_K1_rmse = res_betaZ_K1$rmse, betaZ_K1_bias = res_betaZ_K1$bias,
    betaZ_K2_rmse = res_betaZ_K2$rmse, betaZ_K2_bias = res_betaZ_K2$bias,
    
    PAR_time1 = PAR_time1, PAR_time2 = PAR_time2,
    AAR_time1_k1 = AAR_time1[1], AAR_time1_k2 = AAR_time1[2],
    AAR_time2_k1 = AAR_time2[1], AAR_time2_k2 = AAR_time2[2],
    
    AUC_t1_k1 = auc_t1["t1_k1"], AUC_t1_k2 = auc_t1["t1_k2"],
    AUC_t2_k1 = auc_t2["t2_k1"], AUC_t2_k2 = auc_t2["t2_k2"],
    AUC_t1_k1_mode = auc_t1_mode["t1_k1"], AUC_t1_k2_mode = auc_t1_mode["t1_k2"],
    AUC_t2_k1_mode = auc_t2_mode["t2_k1"], AUC_t2_k2_mode = auc_t2_mode["t2_k2"]
  )
  summary_df <- cbind(summary_df, as.list(flat_alpha), as.list(flat_q))
  list(
    summary_df    = summary_df,
    Q_est_time1   = Q_est_t1, Q_est_time2 = Q_est_t2,
    Q_true_time1  = Q_t1_true, Q_true_time2 = Q_t2_true,
    alpha1_est    = alpha1_pred, alpha2_est = alpha2_pred,
    alpha1_true   = alpha1,      alpha2_true = alpha2,
    diag_df       = diag_df,
    rhat_df       = rhat_df,
    samples       = samples
  )
}


################################################################################
# 4.0 Run one simulation only --------------------------------------------------
################################################################################

#N <- 200
#J <- 6
#r <- 1
#all_runs <- run_one_simulation_rjags_unknownQ(   
#  seed_number = 1000 + 100*N + 10*J + r,
#  N = N, J = J,
#  n.chains = 3, n.adapt = 2000,
#  n.burnin  = 2000, n.iter = 4000, thin = 1  
#)

#results_summary_df <- all_runs$summary_df
#diag_df <- all_runs$diag_df

###############################################################################
# 4. Run Simulations -----------------------------
################################################################################
N_vec        <- c(200, 400, 600)
J_vec        <- 6
replications <- 25

save_every   <- 3
ckpt_file    <- "25_May7_N200400600_J6_RJAGS_unknownQ_1520_final.RDS"

param_grid <- expand.grid(N = N_vec,
                          J = J_vec,
                          r = seq_len(replications),
                          KEEP.OUT.ATTRS = FALSE)

# Load or initialize result container
expected_len <- nrow(param_grid)
if (file.exists(ckpt_file)) {
  all_runs <- readRDS(ckpt_file)
  length(all_runs) <- expected_len   # fill with NULL if needed
} else {
  all_runs <- vector("list", expected_len)
}

for (i in seq_len(expected_len)) {
  if (!is.null(all_runs[[i]])) next
  
  with(param_grid[i, ], {
    cat(sprintf("[N=%d, J=%d] replication %d …\n", N, J, r))
    
    # Repeat until convergence (R-hat < 1.1)
    repeat {
      run_i <- run_one_simulation_rjags_unknownQ(
        seed_number = 1000 + 100*N + 10*J + r,
        N = N, J = J,
        n.chains = 3, n.adapt = 2000,
        n.burnin  = 2000, n.iter = 4000, thin = 1
      )
      max_rhat <- run_i$diag_df$max_rhat_non_alpha
      cat(sprintf("  -> max_rhat_non_alpha = %.3f\n", max_rhat))
      
      if (max_rhat < 1.1) {
        all_runs[[i]] <<- run_i
        break
      } else {
        cat("     R-hat too high, retrying this replication …\n")
      }
    }
  })
  
  if (i %% save_every == 0) {
    saveRDS(all_runs, ckpt_file)
    cat("  progress saved to", ckpt_file, "\n")
  }
}


saveRDS(all_runs, ckpt_file)  # save
cat("✓ All done, full results saved to", ckpt_file, "\n")

#all_runs <- readRDS("25_RJAGS_unknownQ_fullRuns_20250501_1609.RDS")

all_runs <- readRDS("25_May7_N200400600_J6_RJAGS_unknownQ_1520_final.RDS")

# Estract summary_df
summary_list <- lapply(all_runs, \(x) x$summary_df)
results_summary_df <- bind_rows(summary_list)
print(results_summary_df)

diag_all <- do.call(rbind, lapply(all_runs, `[[`, "diag_df"))
print(diag_all)

## 1. AUC ----------
summary_auc <- results_summary_df %>% 
  dplyr::select(N_used, J_used, "AUC_t1_k1", "AUC_t1_k2", 
                "AUC_t2_k1", "AUC_t2_k2") %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_auc)

## 2. AAR ----------
summary_aar <- results_summary_df %>% 
  dplyr::select(N_used, J_used, starts_with("AAR_")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_aar)

## 3. item alpha ----------
summary_gs_A <- results_summary_df %>% 
  dplyr::select(N_used, J_used, g_t1_bias, g_t2_bias, s_t1_bias, s_t2_bias) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_gs_A)

## 4. Q matrix ----------
summary_q <- results_summary_df %>% 
  dplyr::select(N_used, J_used, starts_with("Q_time1_"), starts_with("Q_time2_")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_q)

## 4.1 Metrices of Q matrix ----------
summary_q_t1_metric <- results_summary_df %>% 
  dplyr::select(N_used, J_used, starts_with("q_metrics_t1_")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_q_t1_metric)

summary_q_t2_metric <- results_summary_df %>% 
  dplyr::select(N_used, J_used, starts_with("q_metrics_t2_")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_q_t2_metric)

## 5. bias beta ----------
summary_beta <- results_summary_df %>% 
  dplyr::select(N_used, J_used, beta0_K1_bias, beta0_K2_bias, betaZ_K1_bias, betaZ_K2_bias) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_beta)

## 6. gamma bias ----------
summary_gamma <- results_summary_df %>% 
  dplyr::select(N_used, J_used, gamma01_K1_bias, gamma01_K2_bias, gamma10_K1_bias, gamma10_K2_bias) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_gamma)

## 7. FPR FNR ----------
summary_alpha_FPR <- results_summary_df %>% 
  dplyr::select(N_used, J_used, alpha_metrics_t1_k1_FPR, alpha_metrics_t1_k2_FPR, alpha_metrics_t2_k1_FPR, alpha_metrics_t2_k2_FPR) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_alpha_FPR)

summary_alpha_FNR <- results_summary_df %>% 
  dplyr::select(N_used, J_used, alpha_metrics_t1_k1_FNR, alpha_metrics_t1_k2_FNR, alpha_metrics_t2_k1_FNR, alpha_metrics_t2_k2_FNR) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_alpha_FNR)


# For example:
results_summary_df %>% 
  dplyr::select(N_used, J_used, starts_with("AUC_t1_"), starts_with("AUC_t2_")) %>% 
  group_by(N_used, J_used) %>% 
  reframe(AUC_t1_k1, AUC_t1_k2, AUC_t2_k1, AUC_t2_k2)

results_summary_df %>% 
  dplyr::select(N_used, J_used, starts_with("alpha_metrics_t1_"), starts_with("alpha_metrics_t2_")) %>% 
  group_by(N_used, J_used) %>% 
  reframe(alpha_metrics_t1_k1_FPR,alpha_metrics_t1_k2_FPR, alpha_metrics_t2_k1_FPR,alpha_metrics_t2_k2_FPR)

