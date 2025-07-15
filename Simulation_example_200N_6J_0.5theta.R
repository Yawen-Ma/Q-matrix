
# -------------------------------------------------------
# Step 1. Load required packages
# -------------------------------------------------------
library(cmdstanr)
library(posterior)
library(coda)        # For MCMC diagnostics
library(ggplot2)     # For plotting (if needed)
library(bayesplot)
library(readr)
library(dplyr)
library(MASS)        # for mvrnorm

# -------------------------------------------------------
# Step 2. Define the simulation function
# -------------------------------------------------------

run_one_simulation <- function(seed_number = 123, 
                               model_path = "~/Downloads/Rstudio_my_code/Amplify_games_test/2024new_dataset/2024_newdataset/Q-matrix/20_stan_0.5.stan") {
  
  # --------------------
  # 1) Set the random seed 
  # --------------------
  # This ensures each replication uses a unique seed.
  set.seed(seed_number)
  
  # --------------------
  # 2) Simulate data
  # --------------------
  N <- 200  #263  # Number of students
  J <- 6    # Number of items
  K <- 2    # Number of attributes
  
  # Generate covariates
  Z <- MASS::mvrnorm(N, mu = rep(0, 3), Sigma = diag(3))
  colnames(Z) <- c("Z1","Z2","Z3")
  P <- ncol(Z)
  
  # Generate Q matrices (Time1)
  simulate_Q_time1 <- function(K, J){
    
    Q_time1 <- matrix(0, nrow = J, ncol = K)
    
    #1) identification
    Q_time1[1,] <- c(1,0)
    Q_time1[4,] <- c(0,1)
    
    #2) half-half 1 and 0, since theta = 0.5; each column has one 1
    rows_to_fill <- c(2:3,5:J)
    for (r in rows_to_fill){
      c_chosen = sample(1:2, size = 1)
      Q_time1[r, c_chosen] <- 1
    }
    
    return(Q_time1)
  }
  
  # Generate Q matrices (Time2)
  simulate_Q_time2 <- function(K,J){
    
    Q_time2 <- matrix(0, nrow = J, ncol = K)
    
    #1) identification
    Q_time2[1,] <- c(0,1)
    Q_time2[2,] <- c(1,0)
    
    #2) half-half 1 and 0, since theta = 0.5; each column has one 1
    rows_to_fill <- c(3:J)
    for (r in rows_to_fill){
      c_chosen = sample(1:2, size = 1)
      Q_time2[r, c_chosen] <- 1
    }
    
    return(Q_time2)
  }
  
  Q_time1 <- simulate_Q_time1(K, J)
  Q_time2 <- simulate_Q_time2(K, J)
  
  # Generate true g and s parameters
  gs_true_time1 <- matrix(runif(J*2, 0.05, 0.20), ncol=2)  # slip, guess in [0.05, 0.20]
  gs_true_time2 <- matrix(runif(J*2, 0.05, 0.20), ncol=2)
  
  # Helper functions
  init_probs <- function(z_vec, beta0, betaZ) {
    logit_p <- beta0 + betaZ %*% z_vec
    p <- 1/(1 + exp(-logit_p))
    return(p)
  }
  transition_probs <- function(alpha_prev, z_vec, gamma01, gamma10) {
    p <- numeric(length(alpha_prev))
    z_extended <- c(1, z_vec)  
    for(k in 1:length(alpha_prev)) {
      if(alpha_prev[k] == 0) {
        logit_p <- gamma01[k,] %*% z_extended
        p[k] <- 1/(1 + exp(-logit_p))
      } else {
        logit_p <- gamma10[k,] %*% z_extended
        p[k] <- 1/(1 + exp(-logit_p))
      }
    }
    return(p)
  }
  generate_responses <- function(alpha, Q, gs) {
    N <- nrow(alpha)
    J <- nrow(Q)
    Y <- matrix(NA, nrow=N, ncol=J)
    for(i in 1:N) {
      for(j in 1:J) {
        eta <- prod(alpha[i,]^Q[j,])
        p <- (1-gs[j,2])^eta * gs[j,1]^(1-eta)
        Y[i,j] <- rbinom(1, 1, p)
      }
    }
    return(Y)
  }
  
  #Not good: Generate alpha1
  #beta0 <- rep(0, K)  # intercepts
  #betaZ <- matrix(c(-0.5, 0.3, 0.2,
  #                 0.2, -0.4, 0.1), nrow=K, byrow=TRUE)
  #Better prior
  beta0 <- c(-1, -0.5)
  betaZ <- matrix(c(1.0, 0.8, 0.5,
                    -1.0, -0.5, 0.7), nrow=K, byrow=TRUE)
  
  alpha1 <- matrix(NA, N, K)
  for(i in 1:N) {
    p_init <- init_probs(Z[i,], beta0, betaZ)
    alpha1[i,] <- rbinom(K, 1, p_init)
  }
  
  # Generate transitions
  #gamma01 <- matrix(c(0.2, 0.3, -0.4, 0.1,
  #                    -0.3, 0.4, 0.2, -0.2), nrow=K, ncol=4, byrow=TRUE)
  #gamma10 <- matrix(c(-0.2, -0.3, 0.4, 0.2,
  #                    0.3, -0.4, -0.2, 0.1), nrow=K, ncol=4, byrow=TRUE)
  gamma01 <- matrix(c(0.5, 0.3, -0.4, 0.1,
                      0.8, 0.4, 0.2, -0.2), nrow=K, ncol=4, byrow=TRUE)
  gamma10 <- matrix(c(-0.5, -0.3, 0.4, 0.2,
                      -0.5, -0.4, -0.2, 0.1), nrow=K, ncol=4, byrow=TRUE)
  
  alpha2 <- matrix(NA, N, K)
  for(i in 1:N) {
    p_trans <- transition_probs(alpha1[i,], Z[i,], gamma01, gamma10)
    # alpha2[i,k] depends on alpha1[i,k]
    for(k in 1:K) {
      if(alpha1[i,k] == 0) {
        alpha2[i,k] <- rbinom(1, 1, p_trans[k])
      } else {
        alpha2[i,k] <- rbinom(1, 1, 1 - p_trans[k])
      }
    }
  }
  
  # Generate Y1, Y2
  Y1 <- generate_responses(alpha1, Q_time1, gs_true_time1)
  Y2 <- generate_responses(alpha2, Q_time2, gs_true_time2)
  
  # Prepare stan data
  stan_data <- list(
    N = N,
    J = J,
    K = K,
    Y1 = Y1,
    Y2 = Y2,
    P = P,
    Z = Z
  )
  
  # -----------------------------------------------------
  # 3) Run Stan Model
  # -----------------------------------------------------
  mod <- cmdstan_model(model_path)
  fit <- mod$sample(
    data = stan_data,
    seed = seed_number,     
    chains = 3,
    parallel_chains = 3,
    iter_warmup = 500,
    iter_sampling = 1000,
    adapt_delta = 0.90,
    max_treedepth = 10
  )
  
  # -----------------------------------------------------
  # 4) Extract and compute metrics
  # -----------------------------------------------------
 
  # Diagnostics: Rhat, ESS, divergences, E-BFMI
  # a) Extract draws
  draws_array <- fit$draws()
  
  # b) Summarize draws with posterior::summarise_draws
  diag_summary <- summarise_draws(
    draws_array,
    "mean",        
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  
  # c) Simple check for Rhat threshold
  threshold <- 1.1
  valid_rhat <- diag_summary$rhat[!is.na(diag_summary$rhat)]
  any_high_rhat <- any(valid_rhat > threshold)
  
  # d) Sampler diagnostics
  sampler_diag <- fit$sampler_diagnostics(inc_warmup = FALSE)
  divergent <- sampler_diag[, , "divergent__"]
  n_divergent <- sum(divergent)
  
  # e) E-BFMI
  energy <- sampler_diag[, , "energy__"]
  energy_vec <- as.vector(energy)
  calc_ebfmi_custom <- function(energy) {
    energy <- as.numeric(energy)
    if(length(energy) < 2) return(NA)
    v <- var(energy)
    dE <- diff(energy)
    mean_dE2 <- mean(dE^2)
    ebfmi <- v / mean_dE2
    return(ebfmi)
  }
  ebfmi_value <- calc_ebfmi_custom(energy_vec)
  
  # We'll also store these diagnostic stats in a small data frame
  diag_df <- data.frame(
    seed_used = seed_number,
    max_rhat = max(valid_rhat),
    any_rhat_above_1.1 = any_high_rhat,
    n_divergent = n_divergent,
    ebfmi = ebfmi_value,
    stringsAsFactors = FALSE
  )
  
  # Summaries
  stan_summary <- fit$summary()
  
  # Extract Q_time1, Q_time2 means
  q1_rows <- grep("^Q_time1\\[", stan_summary$variable)
  q2_rows <- grep("^Q_time2\\[", stan_summary$variable)
  
  Q_est_time1_vec <- stan_summary$mean[q1_rows]
  Q_est_time2_vec <- stan_summary$mean[q2_rows]
  Q_est_time1 <- matrix(Q_est_time1_vec, nrow = J, ncol = K, byrow = FALSE)
  Q_est_time2 <- matrix(Q_est_time2_vec, nrow = J, ncol = K, byrow = FALSE)
  
  # Extract item params
  g1_est <- stan_summary$mean[grep("^g1\\[", stan_summary$variable)]
  s1_est <- stan_summary$mean[grep("^s1\\[", stan_summary$variable)]
  g2_est <- stan_summary$mean[grep("^g2\\[", stan_summary$variable)]
  s2_est <- stan_summary$mean[grep("^s2\\[", stan_summary$variable)]
  
  # Extract gamma, beta
  gamma01_est <- stan_summary$mean[grep("^gamma01\\[", stan_summary$variable)]
  gamma10_est <- stan_summary$mean[grep("^gamma10\\[", stan_summary$variable)]
  gamma01_est_matrix <- matrix(gamma01_est, nrow = K, byrow = FALSE)
  gamma10_est_matrix <- matrix(gamma10_est, nrow = K, byrow = FALSE)
  
  beta0_est <- stan_summary$mean[grep("^beta0\\[", stan_summary$variable)]
  betaZ_est <- stan_summary$mean[grep("^betaZ\\[", stan_summary$variable)]
  betaZ_est_matrix <- matrix(betaZ_est, nrow = K, byrow = FALSE)
  
  # Extract alpha posterior means
  mcmc_mat <- as_draws_matrix(fit)
  alpha1_cols <- grep("^alpha1\\[", colnames(mcmc_mat))
  alpha2_cols <- grep("^alpha2\\[", colnames(mcmc_mat))
  alpha1_post_means <- colMeans(mcmc_mat[, alpha1_cols])
  alpha2_post_means <- colMeans(mcmc_mat[, alpha2_cols])
  
  alpha1_est_matrix <- matrix(alpha1_post_means, nrow = N, ncol = K, byrow = FALSE)
  alpha2_est_matrix <- matrix(alpha2_post_means, nrow = N, ncol = K, byrow = FALSE)
  
  alpha1_pred <- ifelse(alpha1_est_matrix > 0.5, 1, 0)
  alpha2_pred <- ifelse(alpha2_est_matrix > 0.5, 1, 0)
  
  # Helper functions to compute PAR, AAR
  compute_PAR <- function(alpha_pred, alpha_true) {
    N <- nrow(alpha_true)
    K <- ncol(alpha_true)
    correct_dims_per_student <- rowSums(alpha_pred == alpha_true)
    pattern_matches <- ifelse(correct_dims_per_student == K, 1, 0)
    par_value <- mean(pattern_matches)
    return(par_value)
  }
  compute_AAR <- function(alpha_pred, alpha_true) {
    N <- nrow(alpha_true)
    K <- ncol(alpha_true)
    AAR_vector <- numeric(K)
    for (k in 1:K) {
      AAR_vector[k] <- mean(alpha_pred[,k] == alpha_true[,k])
    }
    return(AAR_vector)
  }
  
  PAR_time1 <- compute_PAR(alpha1_pred, alpha1)
  AAR_time1 <- compute_AAR(alpha1_pred, alpha1)
  PAR_time2 <- compute_PAR(alpha2_pred, alpha2)
  AAR_time2 <- compute_AAR(alpha2_pred, alpha2)
  
  # Evaluation metrics
  evaluate_binary_matrix <- function(est, true) {
    accuracy <- mean(round(est) == true)
    rmse <- sqrt(mean((est - true)^2))
    list(accuracy = accuracy, rmse = rmse)
  }
  evaluate_continuous_matrix <- function(est, true) {
    rmse <- sqrt(mean((est - true)^2))
    mae  <- mean(abs(est - true))
    bias <- mean(est - true)
    list(rmse = rmse, mae = mae, bias = bias)
  }
  
  # Evaluate Q
  res_q1 <- evaluate_binary_matrix(Q_est_time1, Q_time1)
  res_q2 <- evaluate_binary_matrix(Q_est_time2, Q_time2)
  
  # Evaluate item params
  res_g_t1 <- evaluate_continuous_matrix(g1_est, gs_true_time1[,1])
  res_s_t1 <- evaluate_continuous_matrix(s1_est, gs_true_time1[,2])
  
  res_g_t2 <- evaluate_continuous_matrix(g2_est, gs_true_time2[,1])
  res_s_t2 <- evaluate_continuous_matrix(s2_est, gs_true_time2[,2])
  
  # Evaluate gamma, beta
  res_gamma01_K1 <- evaluate_continuous_matrix(gamma01_est_matrix[1,], gamma01[1,])
  res_gamma01_K2 <- evaluate_continuous_matrix(gamma01_est_matrix[2,], gamma01[2,])
  
  res_gamma10_K1 <- evaluate_continuous_matrix(gamma10_est_matrix[1,], gamma10[1,])
  res_gamma10_K2 <- evaluate_continuous_matrix(gamma10_est_matrix[2,], gamma10[2,])
  
  beta0_est_matrix <- matrix(beta0_est, nrow = K, byrow = FALSE)
  res_beta0_K1 <- evaluate_continuous_matrix(beta0_est_matrix[1], beta0[1])
  res_beta0_K2 <- evaluate_continuous_matrix(beta0_est_matrix[2], beta0[2])
  
  betaZ_est_matrix <- matrix(betaZ_est, nrow = K, byrow = FALSE)
  res_betaZ_K1 <- evaluate_continuous_matrix(betaZ_est_matrix[1,], betaZ[1,])
  res_betaZ_K2 <- evaluate_continuous_matrix(betaZ_est_matrix[2,], betaZ[2,])
  
  # -----------------------------------------------------
  # 5) Return as a result data.frame
  # -----------------------------------------------------
  
  summary_df <- data.frame(
    seed_used = seed_number,
    N_used = N,
    J_used = J,
    
    Q_time1_acc  = res_q1$accuracy,
    Q_time1_rmse = res_q1$rmse,
    Q_time2_acc  = res_q2$accuracy,
    Q_time2_rmse = res_q2$rmse,
    
    g_t1_rmse =  res_g_t1$rmse,
    g_t1_mae  =  res_g_t1$mae,
    g_t1_bias =  res_g_t1$bias,
    
    g_t2_rmse =  res_g_t2$rmse,
    g_t2_mae  =  res_g_t2$mae,
    g_t2_bias =  res_g_t2$bias,
    
    s_t1_rmse =  res_s_t1$rmse,
    s_t1_mae  =  res_s_t1$mae,
    s_t1_bias =  res_s_t1$bias,
    
    s_t2_rmse =  res_s_t2$rmse,
    s_t2_mae  =  res_s_t2$mae,
    s_t2_bias =  res_s_t2$bias,
    
    gamma01_K1_rmse = res_gamma01_K1$rmse,
    gamma01_K1_mae  = res_gamma01_K1$mae,
    gamma01_K1_bias = res_gamma01_K1$bias,
    
    gamma01_K2_rmse = res_gamma01_K2$rmse,
    gamma01_K2_mae  = res_gamma01_K2$mae,
    gamma01_K2_bias = res_gamma01_K2$bias,
    
    gamma10_K1_rmse = res_gamma10_K1$rmse,
    gamma10_K1_mae  = res_gamma10_K1$mae,
    gamma10_K1_bias = res_gamma10_K1$bias,
    
    gamma10_K2_rmse = res_gamma10_K2$rmse,
    gamma10_K2_mae  = res_gamma10_K2$mae,
    gamma10_K2_bias = res_gamma10_K2$bias,
    
    beta0_K1_rmse    = res_beta0_K1$rmse,
    beta0_K1_mae     = res_beta0_K1$mae,
    beta0_K1_bias    = res_beta0_K1$bias,
    
    beta0_K2_rmse    = res_beta0_K2$rmse,
    beta0_K2_mae     = res_beta0_K2$mae,
    beta0_K2_bias    = res_beta0_K2$bias,
    
    betaZ_K1_rmse    = res_betaZ_K1$rmse,
    betaZ_K1_mae     = res_betaZ_K1$mae,
    betaZ_K1_bias    = res_betaZ_K1$bias,
    
    betaZ_K2_rmse    = res_betaZ_K2$rmse,
    betaZ_K2_mae     = res_betaZ_K2$mae,
    betaZ_K2_bias    = res_betaZ_K2$bias,
    
    PAR_time1 = PAR_time1,
    PAR_time2 = PAR_time2,
    
    # Decompose AAR vectors
    AAR_time1_k1 = AAR_time1[1],
    AAR_time1_k2 = AAR_time1[2],
    AAR_time2_k1 = AAR_time2[1],
    AAR_time2_k2 = AAR_time2[2]
  )
  
  return(list(
    summary_df = summary_df,
    Q_est_time1 = Q_est_time1,
    Q_est_time2 = Q_est_time2,
    Q_true_time1 = Q_time1,
    Q_true_time2 = Q_time2,
    
    alpha1_est = alpha1_pred,
    alpha2_est = alpha2_pred,
    alpha1_true = alpha1,   
    alpha2_true = alpha2,
    
    diag_df    = diag_df,
    
    fit = fit
    
  ))
}

# -------------------------------------------------------
# Step 3. Run multiple replications
# -------------------------------------------------------

n_replications <- 50 # CHANGE HERE: if want to have less replications
all_results <- vector("list", n_replications) 

# Help for generating trace plot
all_draws_list <- vector("list", n_replications) 

for (r in 1:n_replications) {
  cat("\n============================\n")
  cat("Starting replication:", r, "\n")
  cat("============================\n")
  
  # We'll setoff the seed by r to create different seeds
  this_seed <- 123 + r
  
  # run_one_simulation now returns a data.frame of results
  sim_output <- run_one_simulation(
    seed_number = this_seed,
    model_path = "~/Downloads/Rstudio_my_code/Amplify_games_test/2024new_dataset/2024_newdataset/Q-matrix/20_stan_0.5.stan"
  )
  
  all_results[[r]] <- sim_output
  
  # store fit to generate trace plot
  fit_object <- sim_output$fit
  
  # select interesting parameters to check trace plot
  draws_selected <- fit_object$draws(variables = c("beta0", "betaZ", "gamma01", "gamma10", "g1", "g2", "s1", "s2", "theta"))
  
  all_draws_list[[r]] <- draws_selected # for r-th replication
}

# Combine all single-row data.frames into one multi-row data.frame
results_summary_df <- do.call(rbind, lapply(all_results, `[[`, "summary_df"))
results_diag_df    <- do.call(rbind, lapply(all_results, `[[`, "diag_df"))

cat("===== Summary Metrics =====\n")
print(results_summary_df)

cat("\n===== Diagnostics =====\n")
print(results_diag_df)

# We can see how many seeds had any Rhat>1.1
cat("\nNumber of replications with any Rhat>1.1:\n")
cat(sum(results_diag_df$any_rhat_above_1.1), " out of ", nrow(results_diag_df), "\n")

cat("\nDone.\n")


# -------------------------------------------------------
# Additionnally: Q-matrix false positive/false negative check
# -------------------------------------------------------

K = 2
J = 6
locked_mask_time1 <- matrix(FALSE, nrow=J, ncol=K)
locked_mask_time2 <- matrix(FALSE, nrow=J, ncol=K)

# Time1
locked_mask_time1[1,1] <- TRUE
locked_mask_time1[1,2] <- TRUE
locked_mask_time1[4,1] <- TRUE
locked_mask_time1[4,2] <- TRUE

# Time2
locked_mask_time2[1,1] <- TRUE
locked_mask_time2[1,2] <- TRUE
locked_mask_time2[2,1] <- TRUE
locked_mask_time2[2,2] <- TRUE

# Define a helper function to compute confusion matrix counts and rates
compute_qmatrix_metrics <- function(Q_est, Q_true, locked_mask, threshold = 0.5) {
  # 1) Ignore identification rows
  idx_free <- which(!locked_mask, arr.ind = TRUE)
  
  # Convert the continuous posterior means to binary predictions (0 or 1)
  Q_est_binary <- ifelse(Q_est > threshold, 1, 0)
  
  # Extract free items, exclude identification elements
  Q_est_free  <- Q_est_binary[idx_free]
  Q_true_free <- Q_true[idx_free]
  
  # Summation over all cells (j,k)
  TP <- sum(Q_est_free == 1 & Q_true_free == 1)
  FP <- sum(Q_est_free == 1 & Q_true_free == 0)
  TN <- sum(Q_est_free == 0 & Q_true_free == 0)
  FN <- sum(Q_est_free == 0 & Q_true_free == 1)
  
  FPR <- if ((FP + TN) == 0) NA else FP / (FP + TN)   # False positive rate
  FNR <- if ((FN + TP) == 0) NA else FN / (FN + TP)   # False negative rate
  accuracy <- (TP + TN) / (TP + FP + TN + FN)
  
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  recall    <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  f1_score  <- if (is.na(precision) || is.na(recall) || (precision + recall == 0)) {
    NA
  } else {
    2 * (precision * recall) / (precision + recall)
  }
  
  return(list(
    TP_count = TP,
    FP_count = FP,
    TN_count = TN,
    FN_count = FN,
    FPR = FPR,
    FNR = FNR,
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1 = f1_score
  ))
}

# Create a data.frame to store the results for each replication
df_q_errors <- data.frame()

# Loop through each simulation replicate in all_results
for (r in seq_along(all_results)) {
  
  # Extract the relevant Q matrices
  Q_est_time1 <- all_results[[r]]$Q_est_time1
  Q_est_time2 <- all_results[[r]]$Q_est_time2
  Q_true_time1 <- all_results[[r]]$Q_true_time1
  Q_true_time2 <- all_results[[r]]$Q_true_time2
  
  # Compute the confusion matrix and metrics for Time 1
  metrics_t1 <- compute_qmatrix_metrics(
    Q_est_time1, 
    Q_true_time1, 
    locked_mask_time1, 
    threshold=0.5
  )
  # Compute the confusion matrix and metrics for Time 2
  metrics_t2 <- compute_qmatrix_metrics(
    Q_est_time2, 
    Q_true_time2, 
    locked_mask_time2, 
    threshold=0.5
  ) 
  
  # Summarize in one row
  row_result <- data.frame(
    replication = r,
    # Time1
    # False positive
    # False negative
    # Overall accuracy
    # Time1:
    TP_count_time1 = metrics_t1$TP_count,
    FP_count_time1 = metrics_t1$FP_count,
    TN_count_time1 = metrics_t1$TN_count,
    FN_count_time1 = metrics_t1$FN_count,
    FPR_time1 = metrics_t1$FPR,
    FNR_time1 = metrics_t1$FNR,
    prec_time1 = metrics_t1$precision,
    recall_time1 = metrics_t1$recall,
    f1_time1 = metrics_t1$f1,
    acc_time1 = metrics_t1$accuracy,
    
    # Time2:
    TP_count_time2 = metrics_t2$TP_count,
    FP_count_time2 = metrics_t2$FP_count,
    TN_count_time2 = metrics_t2$TN_count,
    FN_count_time2 = metrics_t2$FN_count,
    FPR_time2 = metrics_t2$FPR,
    FNR_time2 = metrics_t2$FNR,
    prec_time2 = metrics_t2$precision,
    recall_time2 = metrics_t2$recall,
    f1_time2 = metrics_t2$f1,
    acc_time2 = metrics_t2$accuracy
  )
  
  df_q_errors <- rbind(df_q_errors, row_result)
}

# Results
cat("Q-matrix classification metrics across replications:\n")
print(df_q_errors)

# We can also compute means across replications
df_means <- sapply(df_q_errors[,-1], mean, na.rm=TRUE)   
df_sds   <- sapply(df_q_errors[,-1], sd, na.rm=TRUE)

cat("\nMean of each metric across all replications:\n")
print(df_means)

cat("\nSD of each metric across all replications:\n")
print(df_sds)

# -------------------------------------------------------
# Additionnally: alpha false positive/false negative check
# -------------------------------------------------------

compute_alpha_metrics <- function(alpha_est, alpha_true) {
  # alpha_est, alpha_true: both are N x K 0/1
  alpha_est_vec <- as.vector(alpha_est)  # length N*K
  alpha_true_vec <- as.vector(alpha_true)
  
  TP <- sum(alpha_est_vec == 1 & alpha_true_vec == 1)
  FP <- sum(alpha_est_vec == 1 & alpha_true_vec == 0)
  TN <- sum(alpha_est_vec == 0 & alpha_true_vec == 0)
  FN <- sum(alpha_est_vec == 0 & alpha_true_vec == 1)
  
  # Calculate metrics
  FPR <- if ((FP + TN) == 0) NA else FP/(FP + TN)
  FNR <- if ((FN + TP) == 0) NA else FN/(FN + TP)
  accuracy <- (TP + TN)/(TP + TN + FP + FN)
  precision <- if ((TP + FP)==0) NA else TP/(TP + FP)
  recall <- if ((TP + FN)==0) NA else TP/(TP + FN)
  f1 <- if (is.na(precision) || is.na(recall) || precision+recall==0) {
    NA
  } else {
    2*(precision*recall)/(precision+recall)
  }
  
  return(list(
    alpha_TP_count = TP,
    alpha_FP_count = FP,
    alpha_TN_count = TN,
    alpha_FN_count = FN,
    alpha_FPR = FPR,
    alpha_FNR = FNR,
    alpha_accuracy = accuracy,
    alpha_precision = precision,
    alpha_recall = recall,
    alpha_f1 = f1
  ))
}

# 1) data.frame 
df_alpha_errors <- data.frame()

# 2) compute metrics for each attribute
for (r in seq_along(all_results)) {
  
  alpha1_true_mat <- all_results[[r]]$alpha1_true
  alpha1_est_mat  <- all_results[[r]]$alpha1_est
  alpha2_true_mat <- all_results[[r]]$alpha2_true
  alpha2_est_mat  <- all_results[[r]]$alpha2_est
  
  metrics_a1 <- compute_alpha_metrics(alpha1_est_mat, alpha1_true_mat)
  metrics_a2 <- compute_alpha_metrics(alpha2_est_mat, alpha2_true_mat)
  
  row_alpha <- data.frame(
    replication = r,
    # alpha1 metrics
    alpha1_TP = metrics_a1$alpha_TP_count,
    alpha1_FP = metrics_a1$alpha_FP_count,
    alpha1_TN = metrics_a1$alpha_TN_count,
    alpha1_FN = metrics_a1$alpha_FN_count,
    
    alpha_FPR_a1 = metrics_a1$alpha_FPR,
    alpha_FNR_a1 = metrics_a1$alpha_FNR,
    alpha_prec_a1 = metrics_a1$alpha_precision,
    alpha_recall_a1 = metrics_a1$alpha_recall,
    alpha_f1_a1 = metrics_a1$alpha_f1,
    alpha_acc_a1 = metrics_a1$alpha_accuracy,
    
    # alpha2
    alpha2_TP = metrics_a2$alpha_TP_count,
    alpha2_FP = metrics_a2$alpha_FP_count,
    alpha2_TN = metrics_a2$alpha_TN_count,
    alpha2_FN = metrics_a2$alpha_FN_count,
    
    alpha_FPR_a2 = metrics_a2$alpha_FPR,
    alpha_FNR_a2 = metrics_a2$alpha_FNR,
    alpha_prec_a2 = metrics_a2$alpha_precision,
    alpha_recall_a2 = metrics_a2$alpha_recall,
    alpha_f1_a2 = metrics_a2$alpha_f1,
    alpha_acc_a2 = metrics_a2$alpha_accuracy
  )
  
  df_alpha_errors <- rbind(df_alpha_errors, row_alpha)
}

# 3) Results
cat("Alpha classification metrics across replications:\n")
print(df_alpha_errors)

# 4) Mean
df_means_alpha <- sapply(df_alpha_errors[,-1], mean, na.rm=TRUE)
df_sds_alpha   <- sapply(df_alpha_errors[,-1], sd, na.rm=TRUE)

cat("\nMean of each alpha metric across all replications:\n")
print(df_means_alpha)

cat("\nSD of each alpha metric across all replications:\n")
print(df_sds_alpha)

cat("\nDone.\n")

# -------------------------------------------------------
# Save all results
# -------------------------------------------------------

#save(
#  results_summary_df, results_diag_df, # Bias and diagnostics
#  all_results, # Detailed results
#  df_q_errors, df_means, df_sds, # Q FRP FNP
#  df_alpha_errors, df_means_alpha, df_sds_alpha, # alpha FRP FNP
#  file = "50replics_200_6_0.5.RData" # 50 replications, N = 200, J_t = 6, theta = 0.5
#)

# Do the same for other simulation conditions and name the RData in the same way as replics_N_J_theta.RData 
# For example, 50 repications for N = 400, J = 12, theta = 0.7 should be saved as "50replics_400_12_0.7.RData" 

# -------------------------------------------------------
# Present all results
# -------------------------------------------------------

# Beta parameters ---------------------------------------------------------

library(dplyr)
library(tidyr)
library(xtable)

theta_vec <- c(0.5, 0.7)
N_vec     <- c(200, 400)
J_vec     <- c(6, 12)

results_all <- data.frame(
  beta      = character(),  
  theta     = numeric(),
  N         = numeric(),
  J_t       = numeric(),
  beta0_k1_bias = numeric(),
  beta0_k1_mae  = numeric(),
  beta0_k1_rmse = numeric(),
  beta0_k2_bias = numeric(),
  beta0_k2_mae  = numeric(),
  beta0_k2_rmse = numeric(),
  betaZ_k1_bias = numeric(),
  betaZ_k1_mae  = numeric(),
  betaZ_k1_rmse = numeric(),
  betaZ_k2_bias = numeric(),
  betaZ_k2_mae  = numeric(),
  betaZ_k2_rmse = numeric(),
  stringsAsFactors = FALSE
)

for(th in theta_vec){
  for(n in N_vec){
    for(j in J_vec){
      fname <- sprintf("50replics_%d_%d_%0.1f.RData", n, j, th)
      if(!file.exists(fname)){
        message("Warning: The file not exists", fname)
        next
      }
      load(fname)  
      
      beta0_k1_bias  <- mean(results_summary_df[["beta0_K1_bias"]], na.rm = TRUE)
      beta0_k1_mae   <- mean(results_summary_df[["beta0_K1_mae"]],  na.rm = TRUE)
      beta0_k1_rmse  <- mean(results_summary_df[["beta0_K1_rmse"]], na.rm = TRUE)
      
      beta0_k2_bias  <- mean(results_summary_df[["beta0_K2_bias"]], na.rm = TRUE)
      beta0_k2_mae   <- mean(results_summary_df[["beta0_K2_mae"]],  na.rm = TRUE)
      beta0_k2_rmse  <- mean(results_summary_df[["beta0_K2_rmse"]], na.rm = TRUE)
      
      betaZ_k1_bias  <- mean(results_summary_df[["betaZ_K1_bias"]], na.rm = TRUE)
      betaZ_k1_mae   <- mean(results_summary_df[["betaZ_K1_mae"]],  na.rm = TRUE)
      betaZ_k1_rmse  <- mean(results_summary_df[["betaZ_K1_rmse"]], na.rm = TRUE)
      
      betaZ_k2_bias  <- mean(results_summary_df[["betaZ_K2_bias"]], na.rm = TRUE)
      betaZ_k2_mae   <- mean(results_summary_df[["betaZ_K2_mae"]],  na.rm = TRUE)
      betaZ_k2_rmse  <- mean(results_summary_df[["betaZ_K2_rmse"]], na.rm = TRUE)
      
      results_all <- rbind(
        results_all,
        data.frame(
          beta           = beta_fixed,
          theta          = th,
          N              = n,
          J_t            = j,
          beta0_k1_bias  = beta0_k1_bias,
          beta0_k1_mae   = beta0_k1_mae,
          beta0_k1_rmse  = beta0_k1_rmse,
          beta0_k2_bias  = beta0_k2_bias,
          beta0_k2_mae   = beta0_k2_mae,
          beta0_k2_rmse  = beta0_k2_rmse,
          betaZ_k1_bias  = betaZ_k1_bias,
          betaZ_k1_mae   = betaZ_k1_mae,
          betaZ_k1_rmse  = betaZ_k1_rmse,
          betaZ_k2_bias  = betaZ_k2_bias,
          betaZ_k2_mae   = betaZ_k2_mae,
          betaZ_k2_rmse  = betaZ_k2_rmse,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}
results_long <- results_all %>%
  pivot_longer(
    cols = c(
      beta0_k1_bias, beta0_k1_mae, beta0_k1_rmse,
      beta0_k2_bias, beta0_k2_mae, beta0_k2_rmse,
      betaZ_k1_bias, betaZ_k1_mae, betaZ_k1_rmse,
      betaZ_k2_bias, betaZ_k2_mae, betaZ_k2_rmse
    ),
    names_to    = c("param","metric"),     
    names_pattern = "(beta0_k1|beta0_k2|betaZ_k1|betaZ_k2)_(bias|mae|rmse)",
    values_to   = "value"
  )
results_wide <- results_long %>%
  pivot_wider(
    id_cols     = c(beta, theta, N, J_t, metric),
    names_from  = param,    
    values_from = value     
  ) %>%
  mutate(metric = factor(metric, levels = c("bias","mae","rmse"))) %>%
  arrange(beta, theta, N, J_t, metric)
results_print <-results_wide
results_print <- results_print %>%
  dplyr::select(
    beta, theta, N, J_t, metric,
    beta0_k1, beta0_k2, betaZ_k1, betaZ_k2
  )

results_print





