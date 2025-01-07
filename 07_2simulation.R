set.seed(123)

library(GDINA)
library(edina)
library(LTCDM)
library(MASS)
# Setting parameters -------------------------------------------------

N <- 263 # Number of students
J <- 6 # Number of items
K <- 2 # Number of attributes
S <- 2^K # Number of latent class: 00 01 10 11

# True Q is only used to generate responses, not used for the following analysis
Q_TRUE = array(c(1,1,1,0,1,1,1,1,1,1,1,1,
            0,0,1,1,1,1,1,0,1,1,1,1), dim = c(12,2))

Rep <- 2  # Number of repetitions

run_simulation <- function(N, J, Rep) {
  
# DINA guessing and slipping parameters only used to generate responses
true_gs <- matrix(runif(J*2, 0.2, 0.45), ncol=2)

# Setting covariates ------------------------------------------------------
# Scenario I: three time-constant covariates from a multivariate standard normal distribution
Z <- mvrnorm(N, mu = rep(0, 3), Sigma = diag(3))  # Z is N x 3
colnames(Z) <- c("Z1","Z2","Z3")

# Scenario I assumes the three covarites affected both the initial state and transition probabilities
# Regression coefficients at initial state (reference：00)
# beta0: length 4 for the intercepts of each class (class 1 is reference)
# betaZ: 4 x 3 matrix for the slopes of three covariates for classes {01}, {10}, {11} respectively.
# Note: The first class (00) is reference with beta0[1] and betaZ[1,] effectively 0 in the model.
beta0 <- c(0,  0,   0,  -1)   
betaZ <- matrix(c(-1, 0, -1,
                  0, -0.5, 0.5,
                  -0.5, -1, 0.2,
                  -0.3, 0.1, -0.2), nrow=4, byrow=TRUE)
# Regression coefficients for transitions
# gamma0: intercept terms for the transition probabilities (4x4 matrix)
# gammaZ: slopes for the three covariates (4x4x3 array)
gamma0 <- array(0, dim=c(S,S)) 

# gammaZ extended to 3 dimensions，the 3rd dimension represents the coefficient of the third covariate 
gammaZ <- array(0, dim=c(S,S,3))
# Covariate Z distribution
#Z1 <- rnorm(N, mean=0, sd=1) 
#Z2 <- rnorm(N, mean=0, sd=1)
# e.g. from {00}(prev_state = 1)，transit to {01}  {10} {11}
gamma0[2,1] <- -0.5; gammaZ[2,1,] <- c(0.3, -0.1, 0.2)
gamma0[3,1] <-  0.2; gammaZ[3,1,] <- c(-0.2, 0.4, -0.3)
gamma0[4,1] <-  0.8; gammaZ[4,1,] <- c(0.4, -0.2, 0.1)

# set for the other prev_state
gamma0[2,2] <- -0.3; gammaZ[2,2,] <- c(0.1, 0.2, -0.1)
gamma0[3,2] <-  0.5; gammaZ[3,2,] <- c(0.2, -0.1, 0.3)
gamma0[4,2] <- -0.7; gammaZ[4,2,] <- c(0.3, 0.0, -0.2)

gamma0[2,3] <- -1.0; gammaZ[2,3,] <- c(0.05, -0.05, 0.1)
gamma0[3,3] <-  0.6; gammaZ[3,3,] <- c(-0.1, 0.05, -0.1)
gamma0[4,3] <-  0.3; gammaZ[4,3,] <- c(0.4, -0.3, 0.2)

gamma0[2,4] <- -0.2; gammaZ[2,4,] <- c(0.2, 0.1, 0.05)
gamma0[3,4] <- -0.5; gammaZ[3,4,] <- c(0.3, -0.2, 0.0)
gamma0[4,4] <-  0.7; gammaZ[4,4,] <- c(-0.1, 0.2, -0.3)


# Regression coefficients at initial state (reference：00)
#beta0 <- c(0,  0,  0, -1)   
#betaZ <- c(-1, 0, -1, -0.5)

# Regression coefficients at later state
#gamma0 <- array(0, dim=c(S,S)) # intercept for each row of the transition probability, i.e., {00} transit to one of {00}, {01}, {10}, {11}
#gammaZ <- array(0, dim=c(S,S)) # regression coefficients for latent state s_t except the one that has been designated as the reference category (i.e., s_t = 1)

#gamma0[2,1] <- -0.5; gammaZ[2,1] <- 0.3
#gamma0[3,1] <- 0.2;  gammaZ[3,1] <- -0.2
#gamma0[4,1] <- 0.8;  gammaZ[4,1] <- 0.4

#gamma0[2,2] <- -0.3; gammaZ[2,2] <- 0.1
#gamma0[3,2] <- 0.5;  gammaZ[3,2] <- 0.2
#gamma0[4,2] <- -0.7; gammaZ[4,2] <- 0.3

#gamma0[2,3] <- -1.0; gammaZ[2,3] <- 0.05
#gamma0[3,3] <- 0.6;  gammaZ[3,3] <- -0.1
#gamma0[4,3] <- 0.3;  gammaZ[4,3] <- 0.4

#gamma0[2,4] <- -0.2; gammaZ[2,4] <- 0.2
#gamma0[3,4] <- -0.5; gammaZ[3,4] <- 0.3
#gamma0[4,4] <- 0.7;  gammaZ[4,4] <- -0.1

# Define functions -------------------------------------------------------------------------

## Initial state probabilities
#init_probs <- function(z_vec, beta0, betaZ) {
#  numer2 <- beta0[2] + betaZ[2]*z_vec
#  numer3 <- beta0[3] + betaZ[3]*z_vec
#  numer4 <- beta0[4] + betaZ[4]*z_vec
#  denom <- 1 + exp(numer2) + exp(numer3) + exp(numer4)
#  p1 <- 1 / denom
#  p2 <- exp(numer2)/denom
#  p3 <- exp(numer3)/denom
#  p4 <- exp(numer4)/denom
#  cbind(p1,p2,p3,p4)
#}
init_probs <- function(z_vec, beta0, betaZ) {
  # z_vec is a vector with 3 length
  numer2 <- beta0[2] + sum(betaZ[2, ] * z_vec)
  numer3 <- beta0[3] + sum(betaZ[3, ] * z_vec)
  numer4 <- beta0[4] + sum(betaZ[4, ] * z_vec)
  
  denom <- 1 + exp(numer2) + exp(numer3) + exp(numer4)
  p1 <- 1 / denom
  p2 <- exp(numer2)/denom
  p3 <- exp(numer3)/denom
  p4 <- exp(numer4)/denom
  c(p1,p2,p3,p4)
}

## transition probability
#transition_probs <- function(prev_state, z_vec, gamma0, gammaZ) {
#  numer2 <- gamma0[2,prev_state] + gammaZ[2,prev_state]*z_vec
#  numer3 <- gamma0[3,prev_state] + gammaZ[3,prev_state]*z_vec
# numer4 <- gamma0[4,prev_state] + gammaZ[4,prev_state]*z_vec
#  denom <- 1 + exp(numer2) + exp(numer3) + exp(numer4)
#  p1 <- 1 / denom
#  p2 <- exp(numer2)/denom
#  p3 <- exp(numer3)/denom
#  p4 <- exp(numer4)/denom
#  c(p1,p2,p3,p4)
#}
transition_probs <- function(prev_state, z_vec, gamma0, gammaZ) {

  numer2 <- gamma0[2, prev_state] + sum(gammaZ[2, prev_state, ] * z_vec)
  numer3 <- gamma0[3, prev_state] + sum(gammaZ[3, prev_state, ] * z_vec)
  numer4 <- gamma0[4, prev_state] + sum(gammaZ[4, prev_state, ] * z_vec)
  
  denom <- 1 + exp(numer2) + exp(numer3) + exp(numer4)
  p1 <- 1 / denom
  p2 <- exp(numer2)/denom
  p3 <- exp(numer3)/denom
  p4 <- exp(numer4)/denom
  c(p1,p2,p3,p4)
}


generate_responses <- function(L, Q, gs) {
  N <- length(L)
  J <- nrow(Q)
  Y <- matrix(NA, nrow=N, ncol=J)
  state_by_att <- rbind(c(0,0), # state 1: 00
                        c(0,1), # state 2: 01
                        c(1,0), # state 3: 10
                        c(1,1)) # state 4: 11
  attribute_matrix <- state_by_att[L, ]
  sim <- simGDINA(N = N, 
                  Q = Q, # Q_TRUE[1:6,] and Q_TRUE[7:12,]
                  gs.parm = gs,
                  model = "DINA",
                  attribute = attribute_matrix)
  Y <- GDINA::extract(sim, what = "dat")
 # true.attributes <- GDINA::extract(sim, what = "attribute") # true latent states is the same as L1
}

edina_Q_estimate <- function(data_matrix, burnin = 1000, chain_length = 2000, K = 2){
  model <- edina(data_matrix, k = K, burnin = burnin, chain_length = chain_length)
  model_summary <- summary(model)
  q_matrix_graph <- q_graph(model, binary = TRUE)
  list(
    ModelSummary = model_summary,
    QMatrixGraph = q_matrix_graph,
    edinaModel = model
  )
}

transition_probability <- function(Z1_val, Z2_val, gamma_01, gamma_10) {
  # Define a function from logit to probability
  logit2prob <- function(x) exp(x)/(1+exp(x))
  
  results <- list()
  
  # for each attribute
  for (k in 1:ncol(gamma_01)) {
    # 0->1 transition
    eta_01 <- gamma_01[1, k] + gamma_01[2, k]*Z1_val + gamma_01[3, k]*Z2_val
    p_01 <- logit2prob(eta_01)
    p_00 <- 1 - p_01
    
    # 1->0 transition
    eta_10 <- gamma_10[1, k] + gamma_10[2, k]*Z1_val + gamma_10[3, k]*Z2_val
    p_10 <- logit2prob(eta_10)
    p_11 <- 1 - p_10
    
    # Score the results 
    result <- c("P(0->0)" = p_00, 
                "P(0->1)" = p_01, 
                "P(1->0)" = p_10, 
                "P(1->1)" = p_11)
    
    # Add attribute name
    results[[paste0("Attribute_", k)]] <- result
  }
  
  return(results)
}

# Simulation processes --------------------------------------------------------------

## Initialize --------------------------------------------------------------
# Initialize all parameters to store results across repetitions
PCA_all <- array(NA, dim = c(Rep, 3, 2))   # PCA for year1 and year2
PCV_all <- array(NA, dim = c(Rep, 3, 2))   # PCV for year1 and year2
CEP_all <- list()                          # store CEP for each repetition
true_gs_all <- array(NA, dim = c(Rep, J, 2))
est_gs_all  <- array(NA, dim = c(Rep, J, 2))
difference_tran_matrix_all <- array(NA, dim = c(4,4,Rep))
MAE_all <- numeric(Rep)
MSE_all <- numeric(Rep)


for (r in 1:Rep) {
  # r = 1
  
  true_gs <- matrix(runif(J*2, 0.2, 0.45), ncol=2)
  # store true parameters at this rep
  true_gs_all[r, , ] <- true_gs
  
## Year 1 ------------------------------------------------------------------
  
  # Year 1 initial latent state
  #p_init <- init_probs(Z1, beta0, betaZ)  
  p_init <- t(apply(Z, 1, function(z_i) init_probs(z_i, beta0, betaZ)))
  L1 <- apply(p_init, 1, function(p) sample(1:4, size=1, prob=p))
 # L1 <- apply(p_init, 1, function(p) sample(1:4, replace=TRUE, prob=p))
  state_by_att <- rbind(c(0,0), # state 1: 00
                        c(0,1), # state 2: 01
                        c(1,0), # state 3: 10
                        c(1,1)) # state 4: 11
  # Make L1 into a matrix
  true.L1 <- state_by_att[L1, ] 
  # Base on L1 to generate Year 1 responses Y1 
  Y1 <- generate_responses(L1, Q_TRUE[1:6,], true_gs)
  
  # Based on Y1, estimate Q matrix
  q_est_result_Y1 <- edina_Q_estimate(Y1, burnin = 1000, chain_length = 2000, K = 2) 
  q1 <- extract_q_matrix(q_est_result_Y1$ModelSummary)
  est_y1 <- GDINA(Y1, q1, model = "DINA", verbose = 0) # fit the model and get the estimations
  # different estimates for student attribute profile
  MLE <- personparm(est_y1, what = "MLE")
  MAP <- personparm(est_y1, what = "MAP")
  EAP <- personparm(est_y1, what = "EAP")
  # PCA (Proportion Correct Assignment)：The proportion of the estimates vectors are the same as the true
  PCA_y1 <- c(ClassRate(EAP, true.L1)$PCA,
              ClassRate(MLE[, 1:K], true.L1)$PCA,
              ClassRate(EAP[, 1:K], true.L1)$PCA)
  # PCV (Proportion Correct Value)
  PCV_y1 <- c(ClassRate(EAP, true.L1)$PCV[K],
              ClassRate(MLE[, 1:K], true.L1)$PCV[K],
              ClassRate(EAP[, 1:K], true.L1)$PCV[K])
  
  # Estimater parameters: guess and slip
  est_gs <- coef(est_y1, what = "gs") 
  est_gs_all[r, , ] <- est_gs
  

## Year 2 ------------------------------------------------------------------
  L2 <- numeric(N)
  for(i in 1:N) {
   # p_trans <- transition_probs(L1[i], Z2[i], gamma0, gammaZ)
   # L2[i] <- sample(1:4, replace=TRUE, prob=p_trans)
    p_trans <- transition_probs(L1[i], Z[i,], gamma0, gammaZ)
    L2[i] <- sample(1:4, size=1, prob=p_trans)
   }
  state_by_att <- rbind(c(0,0), # state 1: 00
                        c(0,1), # state 2: 01
                        c(1,0), # state 3: 10
                        c(1,1)) # state 4: 11
  # Make L1 into a matrix
  true.L2 <- state_by_att[L2, ] 
  
  Y2 <- generate_responses(L2, Q_TRUE[7:12,], true_gs)
  
  q_est_result_Y2 <- edina_Q_estimate(Y2, burnin = 1000, chain_length = 2000, K = 2)
  
  q2 <- extract_q_matrix(q_est_result_Y2$ModelSummary)
  est_y2 <- GDINA(Y2, q2, model = "DINA", verbose = 0) # fit the model and get the estimations
  # different estimates for student attribute profile
  MLE <- personparm(est_y2, what = "MLE")
  MAP <- personparm(est_y2, what = "MAP")
  EAP <- personparm(est_y2, what = "EAP")
  # PCA (Proportion Correct Assignment)：The proportion of the estimates vectors are the same as the true
  PCA_y2 <- c(ClassRate(EAP, true.L2)$PCA,
              ClassRate(MLE[, 1:K], true.L2)$PCA,
              ClassRate(EAP[, 1:K], true.L2)$PCA)
  PCV_y2 <- c(ClassRate(EAP, true.L2)$PCV[K],
              ClassRate(MLE[, 1:K], true.L2)$PCV[K],
              ClassRate(EAP[, 1:K], true.L2)$PCV[K])
  
  # Estimater parameters: guess and slip
  est_gs <- coef(est_y1, what = "gs") 
  est_gs_all[r, , ] <- est_gs
  
## Store PCA and PCV for both years in PCA_all and PCV_all --------
  # dimension: Rep x 3 x 2, where 2 is year index: 1 for year1, 2 for year2
  PCA_all[r, , 1] <- PCA_y1
  PCV_all[r, , 1] <- PCV_y1
  PCA_all[r, , 2] <- PCA_y2
  PCV_all[r, , 2] <- PCV_y2
  
  
  # True transition matrix
  state_counts <- table(L1)
  
  true_transition_count <- matrix(0, nrow=4, ncol=4)
  
  for (i in 1:N) {
    p_trans <- transition_probs(L1[i], Z2[i], gamma0, gammaZ)
    # p_trans is the true transition probability from L1[i] to the other 4 latent states at next time point
    true_transition_count[L1[i], ] <- true_transition_count[L1[i], ] + p_trans
  }
  
  # average the row, considering the row sum up to 1
  for (s in 1:4) {
    if (state_counts[s] > 0) {
      true_transition_count[s, ] <- true_transition_count[s, ] / state_counts[s]
    }
  }
  
  
  true_transition_matrix <- true_transition_count
  rownames(true_transition_matrix) <- c("{00}", "{01}", "{10}", "{11}")
  colnames(true_transition_matrix) <- c("{00}", "{01}", "{10}", "{11}")
  
  # considering the row sum up to 1
  for (s in 1:4) {
    row_sum <- sum(true_transition_matrix[s, ])
    if(row_sum != 1) {
      true_transition_matrix[s, ] <- true_transition_matrix[s, ] / row_sum
    }
  }
  
  
  # Estimate transition matrix
  # Step 1: Put the G-DINA model objects at pre- and post-tests in a list
  fit.object = list()
  fit.object[[1]] <- est_y1
  fit.object[[2]] <- est_y2
  Q_list = list(Q1 = Q_TRUE[1:6,], Q2 = Q_TRUE[7:12,])
  
  # Step 2: Compute classification error probabilities (CEP)
  # The classification error probabilities (CEP) can be obtained in this data example
  cep = cep_t(fit.object = fit.object, t = 2, K = 2, N = 263)
  
##  Store CEP ------
  CEP_all[[r]] <- cep  
  
  # The CEP matrices of the attributes
  cep$cep.matrix
 
  cep$cep.matrix[[1]] # Attribute 1
  cep$cep.matrix[[2]] # Attribute 2
  
  # Step 3: Covariates
#  Z = cbind(Z1, Z2)
  
  z_t1 = cbind(1, Z)
 # z_t2 = cbind(1, Z[,"Z1"], Z[,"Z2"], apply(Z,1,prod)) # # Covariates at time 2
  z_t2 = cbind(1, Z) # # Covariates at time 2
  
  # Set appropriate initial values of the coefficients
  # Initial values for the initial state's regression coefficients
  beta_in = matrix(0, ncol(z_t1), K) 
  # ga01_in: transition from 0 to 1
  #ga01_in = cbind(
  #  c(0,  -0.5, 0, 0, 0, 1),  # Attribute 1
  #  c(-0.1,  0.4, -2.3, -0.9, 0, 0.3)   # Attribute 2
  #)
  ga01_in = cbind(
      c(0,  -0.5, 0),  # Attribute 1
      c(-0.1,  0.4, 0)   # Attribute 2
    )
  ga01_in = cbind(
    c(0.2,  -0.2, 0.3),  # Attribute 1
    c(-0.2,  0.2, 0.1)   # Attribute 2
  )
  # ga10_in: transition from 1 to 0
  # 1->0 good covariate decrease the probability of transition(negative)
  # 1->0 bad covariate increase the probability of transition(positive)）
 # ga10_in = cbind(
  #  c(0.2, -0.2,  -0.3, 0.3,  -0.2, 0.1),  #  Attribute 1
  #  c(-0.2, 0.2,  -0.9, -0.1,  0.3, -0.5)  # Attribute 2
  #)
  # Note that this estimation can take several minutes.
  step3.output <- step3.est(cep = cep, z_t1 = z_t1, z_t2 = z_t2, K = K, t = t, beta_in, ga01_in, ga10_in)
  step3.output
  beta = step3.output$beta
  gamma_01 = step3.output$gamma_01
  gamma_10 = step3.output$gamma_10
  
  # for each user, calculate transition probability
  user_probs_list <- apply(Z, 1, function(row) {
    Z1_val <- as.numeric(row["Z1"])
    Z2_val <- as.numeric(row["Z2"])
    transition_probability(Z1_val, Z2_val, gamma_01, gamma_10)
  })
  
  # attract attribute 1 matrix, dim = 4 x 263
  attr1_matrix <- sapply(user_probs_list, function(x) x[["Attribute_1"]])
  # attract attribute 2 matrix, dim = 4 x 263
  attr2_matrix <- sapply(user_probs_list, function(x) x[["Attribute_2"]])
  
  # Calculate the mean value
  group_means_attr1 <- rowMeans(attr1_matrix)
  group_means_attr2 <- rowMeans(attr2_matrix)
  
  # Rearrange the transition matrix
  p00_1 <- group_means_attr1["P(0->0)"]
  p01_1 <- group_means_attr1["P(0->1)"]
  p10_1 <- group_means_attr1["P(1->0)"]
  p11_1 <- group_means_attr1["P(1->1)"]
  
  p00_2 <- group_means_attr2["P(0->0)"]
  p01_2 <- group_means_attr2["P(0->1)"]
  p10_2 <- group_means_attr2["P(1->0)"]
  p11_2 <- group_means_attr2["P(1->1)"]
  
  joint_matrix <- matrix(nrow=4, ncol=4, byrow=TRUE, data=c(
    # {00}->...
    p00_1*p00_2, p00_1*p01_2, p01_1*p00_2, p01_1*p01_2,
    # {01}->...
    p00_1*p10_2, p00_1*p11_2, p01_1*p10_2, p01_1*p11_2,
    # {10}->...
    p10_1*p00_2, p10_1*p01_2, p11_1*p00_2, p11_1*p01_2,
    # {11}->...
    p10_1*p10_2, p10_1*p11_2, p11_1*p10_2, p11_1*p11_2
  ))
  
  rownames(joint_matrix) <- c("{00}", "{01}", "{10}", "{11}")
  colnames(joint_matrix) <- c("{00}", "{01}", "{10}", "{11}")
  
  # Make sure the sum of each row is 1
  for (s in 1:4) {
    row_sum <- sum(joint_matrix[s, ])
    if (row_sum != 1) {
      joint_matrix[s, ] <- joint_matrix[s, ] / row_sum
    }
  }
  
  difference_tran_matrix = joint_matrix - true_transition_matrix
  difference <- difference_tran_matrix 
  MAE <- mean(abs(difference))
  MSE <- mean((difference)^2)
  
## Store difference_tran_matrix, MAE, MSE for each repetition ------
  difference_tran_matrix_all[,,r] <- difference_tran_matrix
  MAE_all[r] <- MAE
  MSE_all[r] <- MSE
  
  cat("Rep:", r, "\n")
  cat("Mean Absolute Error:", MAE, "\n")
  cat("Mean Squared Error:", MSE, "\n")
  
 # e.g. cat("Mean Absolute Error:", MAE, "\n") # 0.197: This indicates that the estimated transition probabilities deviate from the true values by an average of approximately 20 percentage.
 # e.g. cat("Mean Squared Error:", MSE, "\n") # An MSE of 0.081 implies that the average squared error is 0.08, corresponding to square root of MSE of approximately 0.285.
}



# You can now summarize results:
mean_MAE <- mean(MAE_all)
mean_MSE <- mean(MSE_all)

cat("Average MAE across reps:", mean_MAE, "\n")
cat("Average MSE across reps:", mean_MSE, "\n")

# Similarly, you can look at PCA and PCV averages:
mean_PCA_year1 <- apply(PCA_all[,,1], 2, mean, na.rm=TRUE)
mean_PCA_year2 <- apply(PCA_all[,,2], 2, mean, na.rm=TRUE)
mean_PCV_year1 <- apply(PCV_all[,,1], 2, mean, na.rm=TRUE)
mean_PCV_year2 <- apply(PCV_all[,,2], 2, mean, na.rm=TRUE)

cat("Mean PCA Year 1:", mean_PCA_year1, "\n")
cat("Mean PCA Year 2:", mean_PCA_year2, "\n")
cat("Mean PCV Year 1:", mean_PCV_year1, "\n")
cat("Mean PCV Year 2:", mean_PCV_year2, "\n")

guess_bias <- numeric(J)
guess_rmse <- numeric(J)
slip_bias <- numeric(J)
slip_rmse <- numeric(J)

for (j in 1:J) {
  true_guess_j <- true_gs_all[, j, 1]
  est_guess_j  <- est_gs_all[, j, 1]
  
  true_slip_j <- true_gs_all[, j, 2]
  est_slip_j  <- est_gs_all[, j, 2]
  
  guess_bias[j] <- mean(est_guess_j - true_guess_j, na.rm = TRUE)
  guess_rmse[j] <- sqrt(mean((est_guess_j - true_guess_j)^2, na.rm = TRUE))
  
  slip_bias[j] <- mean(est_slip_j - true_slip_j, na.rm = TRUE)
  slip_rmse[j] <- sqrt(mean((est_slip_j - true_slip_j)^2, na.rm = TRUE))
}

cat("Guess Bias per item:", guess_bias, "\n")
cat("Guess RMSE per item:", guess_rmse, "\n")
cat("Slip Bias per item:", slip_bias, "\n")
cat("Slip RMSE per item:", slip_rmse, "\n")


# add timestamp
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

# save results
saveRDS(list(
  PCA_all = PCA_all,
  PCV_all = PCV_all,
  CEP_all = CEP_all,
  MAE_all = MAE_all,
  MSE_all = MSE_all
), file = paste0("simulation_results_case1_", N, "N_", J, "J_", Rep, "Rep_", timestamp, ".rds"))

cat("Simulation completed and results saved.\n")

}

## plot ------
library(ggplot2)
library(reshape2) 

item_index <- 1:J
guess_data <- data.frame(
  Item = factor(item_index),
  Bias = guess_bias,
  RMSE = guess_rmse
)
slip_data <- data.frame(
  Item = factor(item_index),
  Bias = slip_bias,
  RMSE = slip_rmse
)

guess_melted <- melt(guess_data, id.vars = "Item", variable.name = "Metric", value.name = "Value")
slip_melted <- melt(slip_data, id.vars = "Item", variable.name = "Metric", value.name = "Value")

# Guess Bias and RMSE
guess_plot <- ggplot(guess_melted, aes(x = Item, y = Value, fill = Metric)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title = "Guess Parameter Bias and RMSE by Item",
       x = "Item",
       y = "Value") +
  scale_fill_discrete(name = "Guess Parameter",
                      labels = c("Bias", "RMSE"))

# Slip parameters
slip_plot <- ggplot(slip_melted, aes(x = Item, y = Value, fill = Metric)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title = "Slip Parameter Bias and RMSE by Item",
       x = "Item",
       y = "Value") +
  scale_fill_discrete(name = "Slip Parameter",
                      labels = c("Bias", "RMSE"))

