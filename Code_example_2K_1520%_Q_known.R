
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
library(pROC) # ROC AUC

# Stan model input here 

mod <- cmdstan_model("24_stan_1520%_Q_known.stan")   # Stan file name

# -------------------------------------------------------
# Step 2. Define the simulation function
# -------------------------------------------------------
# ------------------ Q-matrix  ------------------
make_Q_pair_fixed <- function(J){
  stopifnot(J %in% c(6, 18, 30))
  K <- 2
  
  if(J == 6){
    Q1 <- rbind(
      diag(1, K), diag(1, K),    
      c(1,0),
      c(1,1) 
    )
    Q2 <- rbind(
      c(0,1), c(1,0), 
      c(0,1), c(1,0), 
      c(0,1), c(1,1) 
    )
  }
  
  if(J == 18){
    Q1 <- rbind(
      matrix(rep(diag(1, K), 3), ncol = K, byrow = TRUE), 
      matrix(rep(diag(1, K), 3), ncol = K, byrow = TRUE),
      c(1,0), c(0,1),
      matrix(rep(c(1,1), 4), 4, K, byrow = TRUE)         
    )
    Q2 <- rbind(
      c(0,1), c(1,0),
      matrix(rep(c(0,1), 6), ncol = K, byrow = TRUE),
      matrix(rep(c(1,0), 6), ncol = K, byrow = TRUE),
      matrix(rep(c(1,1), 4), ncol = K, byrow = TRUE)
    )
  }
  
  if(J == 30){
    Q1 <- rbind(
      matrix(rep(diag(1, K), 6), ncol = K, byrow = TRUE), 
      matrix(rep(diag(1, K), 6), ncol = K, byrow = TRUE),
      matrix(rep(c(1,1), 6), 6, K, byrow = TRUE)
    )
    Q2 <- rbind(
      c(0,1), c(1,0),
      matrix(rep(c(0,1), 11), ncol = K, byrow = TRUE),
      matrix(rep(c(1,0), 11), ncol = K, byrow = TRUE),
      matrix(rep(c(1,1), 6), ncol = K, byrow = TRUE)
    )
  }
  
  stopifnot(nrow(Q1) == J, nrow(Q2) == J)
  list(Q_time1 = Q1, Q_time2 = Q2)
}


# ----------  Q / alpha ---------- 
compute_confusion <- function(est, true){
  v_est  <- as.integer(est>0.5); v_true <- as.integer(true)
  TP <- sum(v_est & v_true)
  FP <- sum(v_est & !v_true)
  TN <- sum(!v_est & !v_true)
  FN <- sum(!v_est &  v_true)
  data.frame(
    TP = TP, FP = FP, TN = TN, FN = FN,
    FPR = ifelse((FP+TN)==0, NA, FP/(FP+TN)),
    FNR = ifelse((FN+TP)==0, NA, FN/(FN+TP)),
    TPR = ifelse((TP+FN)==0, NA, TP/(TP+FN)),  # True Positive Rate
    TNR = ifelse((TN+FP)==0, NA, TN/(TN+FP)),   # True Negative Rate
    ACC = (TP+TN)/(TP+FP+TN+FN)
  )
}

run_one_simulation <- function(seed_number = 123, 
                               N          = 200,
                               J          = 6,
                               model_path = "~/Downloads/Rstudio_my_code/2024_newdataset/Q-matrix/24_stan_1520%_Q_known.stan" # Change into your address of the stan file
) {
  
  # --------------------
  # 1) Set the random seed 
  # --------------------
  # This ensures each replication uses a unique seed.
  set.seed(seed_number)
  
  # --------------------
  # 2) Simulate data
  K <- 2    # Number of attributes
  
  # Generate covariates
  Z <- MASS::mvrnorm(N, mu = rep(0, 3), Sigma = diag(3))
  colnames(Z) <- c("Z1","Z2","Z3")
  P <- ncol(Z)
  
  Qs <- make_Q_pair_fixed(J)
  Q_time1 <- Qs$Q_time1
  Q_time2 <- Qs$Q_time2
  
  set.seed(123)
  
  gs_true_time1 <- matrix(runif(J*2, 0.05, 0.15), ncol=2)
  gs_true_time2 <- matrix(runif(J*2, 0.05, 0.15), ncol=2)
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

  set.seed(123) 
  beta0 <- rnorm(K, mean = 0, sd = 0.6)  # beta0 ~ N(0, 2)
  betaZ <- matrix(rnorm(K * 3, mean =  0, sd = 0.3), nrow = K, ncol = 3)
 
  # time point 1
  alpha1 <- matrix(NA, N, K)
  
  for(i in 1:N) {
    p_init <- init_probs(Z[i,], beta0, betaZ)
    alpha1[i,] <- rbinom(K, 1, p_init)
  }
  
  alpha1 %>% as.data.frame(alpha1) %>% count(V1, V2)
  set.seed(123) 

  gamma01 <- matrix(NA, nrow = K, ncol = P+1)  
  gamma10 <- matrix(NA, nrow = K, ncol = P+1)
  
  gamma01[1,1] <- rnorm(1, mean = -1.66, sd = 0.3)
  gamma10[1,1] <- rnorm(1, mean = -4.87, sd = 1)
  gamma01[2,1] <- rnorm(1, mean = -2.0, sd = 0.3)
  gamma10[2,1] <- rnorm(1, mean = -2.2, sd = 0.3)
  
  for(k in 1:2) {
    for(pp in 2:(P+1)) {
      gamma01[k,pp] <- rnorm(1, mean = 0, sd = 0.3)
      gamma10[k,pp] <- rnorm(1, mean = 0, sd = 0.3)
    }
  }
  
  # time point 2
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
  alpha2 %>% as.data.frame(alpha1) %>% count(V1, V2)
  
  # Generate Y1, Y2
  Y1 <- generate_responses(alpha1, Q_time1, gs_true_time1)
  Y2 <- generate_responses(alpha2, Q_time2, gs_true_time2)
  
  stan_data <- list(
    N = N,
    J = J,
    Y1 = Y1,
    Y2 = Y2,
    P = P,
    Z = Z,
    Q1 = Q_time1, 
    Q2 = Q_time2
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
    max_treedepth = 10,
    save_cmdstan_config = TRUE)
  
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
  alpha1_cols <- grep("^prob_attr1\\[", colnames(mcmc_mat))
  alpha2_cols <- grep("^prob_attr2\\[", colnames(mcmc_mat))
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
  
  AAR_time1
  AAR_time2
 
  # ------------------------------------------------------------
  ## 1) helpful function for AUC
  ## ------------------------------------------------------------
  get_auc_vec <- function(prob_mat, truth_mat, prefix = "AUC") {
    K  <- ncol(truth_mat)
    au <- numeric(K)
    for (k in seq_len(K)) {
      roc_k <- roc(response  = truth_mat[, k],     # 0/1 
                   predictor = prob_mat[, k],      
                   quiet     = TRUE,
                   direction = "<")                
      au[k] <- as.numeric(auc(roc_k))
    }
    names(au) <- sprintf("%s_k%d", prefix, seq_len(K))
    au
  }
  
  ## ------------------------------------------------------------
  ## 2) AUC for time 1 and time 2
  ## ------------------------------------------------------------
  AUC_time1 <- get_auc_vec(alpha1_est_matrix, alpha1, prefix = "t1")
  AUC_time2 <- get_auc_vec(alpha2_est_matrix, alpha2, prefix = "t2")
  
  AUC_time1
  AUC_time2
  AUC_t1_k1 = AUC_time1["t1_k1"]
  AUC_t1_k2 = AUC_time1["t1_k2"]
  AUC_t2_k1 = AUC_time2["t2_k1"]
  AUC_t2_k2 = AUC_time2["t2_k2"]
  
   # Evaluation metrics
  evaluate_binary_matrix <- function(est, true) {
    accuracy <- mean(round(est) == true)
    rmse <- sqrt(mean((est - true)^2))
    list(accuracy = accuracy, rmse = rmse)
  }
  evaluate_continuous_matrix <- function(est, true) {
    rmse <- sqrt(mean((est - true)^2))
    bias <- mean(est - true)
    list(rmse = rmse, bias = bias)
  }
  
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
  
  ## ---------- (4) confusion-matrix ----------
  alpha_metrics_t1_k1 <- compute_confusion(alpha1_pred[,1], alpha1[,1])
  alpha_metrics_t1_k2 <- compute_confusion(alpha1_pred[,2], alpha1[,2])
  alpha_metrics_t2_k1 <- compute_confusion(alpha2_pred[,1], alpha2[,1])
  alpha_metrics_t2_k2 <- compute_confusion(alpha2_pred[,2], alpha2[,2])
  
  alpha_list <- list(
    t1_k1 = alpha_metrics_t1_k1,
    t1_k2 = alpha_metrics_t1_k2,
    t2_k1 = alpha_metrics_t2_k1,
    t2_k2 = alpha_metrics_t2_k2
  )
  flatten_confusion <- function(df, prefix){
    metrics <- c("FPR","FNR","TPR","TNR","ACC")
    vals    <- as.numeric(df[1, metrics])
    setNames(vals, paste0("alpha_metrics_", prefix, "_", metrics))
  }
  flat_alpha <- unlist(
    lapply(names(alpha_list),
           \(nm) flatten_confusion(alpha_list[[nm]], nm))
  )
  
  # -----------------------------------------------------
  # 5) Return as a result data.frame
  # -----------------------------------------------------
  
  summary_df <- data.frame(
    seed_used = seed_number,
    N_used = N,
    J_used = J,
    
    g_t1_rmse =  res_g_t1$rmse,
    g_t1_bias =  res_g_t1$bias,
    
    g_t2_rmse =  res_g_t2$rmse,
    g_t2_bias =  res_g_t2$bias,
    
    s_t1_rmse =  res_s_t1$rmse,
    s_t1_bias =  res_s_t1$bias,
    
    s_t2_rmse =  res_s_t2$rmse,
    s_t2_bias =  res_s_t2$bias,
    
    gamma01_K1_rmse = res_gamma01_K1$rmse,
    gamma01_K1_bias = res_gamma01_K1$bias,
    
    gamma01_K2_rmse = res_gamma01_K2$rmse,
    gamma01_K2_bias = res_gamma01_K2$bias,
    
    gamma10_K1_rmse = res_gamma10_K1$rmse,
    gamma10_K1_bias = res_gamma10_K1$bias,
    
    gamma10_K2_rmse = res_gamma10_K2$rmse,
    gamma10_K2_bias = res_gamma10_K2$bias,
    
    beta0_K1_rmse    = res_beta0_K1$rmse,
    beta0_K1_bias    = res_beta0_K1$bias,
    
    beta0_K2_rmse    = res_beta0_K2$rmse,
    beta0_K2_bias    = res_beta0_K2$bias,
    
    betaZ_K1_rmse    = res_betaZ_K1$rmse,
    betaZ_K1_bias    = res_betaZ_K1$bias,
    
    betaZ_K2_rmse    = res_betaZ_K2$rmse,
    betaZ_K2_bias    = res_betaZ_K2$bias,
    
    PAR_time1 = PAR_time1,
    PAR_time2 = PAR_time2,
    
    # Decompose AAR vectors
    AAR_time1_k1 = AAR_time1[1],
    AAR_time1_k2 = AAR_time1[2],
    AAR_time2_k1 = AAR_time2[1],
    AAR_time2_k2 = AAR_time2[2],
    
    AUC_t1_k1 = AUC_t1_k1,
    AUC_t1_k2 = AUC_t1_k2,
    AUC_t2_k1 = AUC_t2_k1,
    AUC_t2_k2 = AUC_t2_k2
    )
  
  summary_df <- cbind(summary_df, as.list(flat_alpha))
  
  return(list(
    summary_df = summary_df,
    
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
N_vec         <- c(200, 600)#, 1000)
J_vec         <- c(6, 18)#, 30)
replications  <- 10

all_runs <- list(); row <- 1

for(N in N_vec){
  for(J in J_vec){
    for(r in seq_len(replications)){
      cat(sprintf("[N=%d, J=%d] replication %d ...\n", N, J, r))
      res <- run_one_simulation(N = N,
                                J = J,
                                seed = 1000 + 100*N + 10*J + r)
      all_runs[[row]] <- res$summary
      row <- row + 1
    }
  }
}

results_summary_df_1520 <- do.call(rbind, all_runs)

print(results_summary_df_1520)

out_file <- sprintf(
  "23_2K_1520_Q_known_pct_27Apr_SimSummary_%s.RDS",          
  format(Sys.time(), "%Y%m%d_%H%M")
)

saveRDS(results_summary_df_1520, out_file)
cat("✓ summary saved to", out_file, "\n")

results_summary_df_1520 <- readRDS("23_2K_1520pct_27Apr_SimSummary_20250427_0726.RDS")

## 1. all metrics mean/sd ----------
summary_stats <- results_summary_df_1520 %>% 
  group_by(N_used, J_used) %>%                    
  summarise(across(
    .cols  = -seed_used,                     
    .fns   = list(mean = mean, sd = sd),
    .names = "{.col}_{.fn}"
  ), .groups = "drop")

print(summary_stats)

## 2. PAR/AAR ----------
summary_par_aar <- results_summary_df_1520 %>% 
  dplyr::select(N_used, J_used, starts_with("PAR_"), starts_with("AAR_")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_par_aar)

## 3. bias alpha ----------
summary_alpha <- results_summary_df_1520 %>% 
  dplyr::select(N_used, J_used, starts_with("bias")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean, sd = sd)),
            .groups = "drop")

print(summary_alpha)

## 4. item alpha ----------
summary_gs_A <- results_summary_df_1520 %>% 
  dplyr::select(N_used, J_used, starts_with("g_"), starts_with("s_")) %>% 
  group_by(N_used, J_used) %>% 
  summarise(across(everything(), list(mean = mean)),
            .groups = "drop")

print(summary_gs_A)

