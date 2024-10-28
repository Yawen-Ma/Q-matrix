# Estimate K and Q --------------------------------------------------------

# We have different packages available to estimate the K (the number of attribute)
# and Q-matrix 

# Introduction to packages --------------------------------------------------------

# {edina}: developed in 2020 by James Joseph Balamuta, Steven Andrew Culpepper, Jeffrey A. Douglas https://cran.r-project.org/web/packages/edina/edina.pdf
# relate with paper in 2018 "Bayesian estimation of the DINA Q-matrix" (Y. Chen, S.A. Culpepper, Y. Chen, J.A. Douglas) see https://link.springer.com/article/10.1007/s11336-017-9579-4

# {edina} --------------------------------------------------------

# Set a seed for reproducibility of random operations
set.seed(123)

## {edina} K estimate --------------------------------------------------------

edina_K_estimate <- function(cluster_subset){
  # Convert data to matrix and remove the first column if it's an ID or non-numeric
  data_matrix <- as.matrix(cluster_subset[,-1])
  
  # Running the EDINA model
  model <- auto_edina(data_matrix, burnin = 1000, chain_length = 2000)
  
  # Generate summary of the model
  model_summary <- summary(model)
  
  # Identify the best model configuration
  best_model_result <- best_model(model)
  
  # Return all results in a list
  return(list(
    Summary = model_summary,
    BestModel = best_model_result
  ))
}

# For example, cluster_subset is C5_subset
C5_edina_K <- edina_K_estimate(C5_subset)

print(C5_edina_K)

# The EDINA model for data with K = 3/4 

## {edina} Q estimate --------------------------------------------------------

edina_Q_estimate <- function(cluster_subset, burnin = 1000, chain_length = 2000){
  # Preparing the matrix for the EDINA model, removing the first column if it's not relevant
  data_matrix <- as.matrix(cluster_subset[, -1])
  
  # Run the EDINA model, set K = 4 manually
  model <- edina(data_matrix, k = 4, burnin = burnin, chain_length = chain_length)
  
  # Generate the summary of the model
  model_summary <- summary(model)
  
  # Visual representation of the Q matrix
  q_matrix_graph <- q_graph(model, binary = TRUE)
  
  # Returning a list of all the results for further analysis or reporting
  return(list(
    ModelSummary = model_summary,
    QMatrixGraph = q_matrix_graph
  ))
  
}

C5_edina_Q <- edina_Q_estimate(C5_subset)

print(C5_edina_Q)

# To compare with results from cdmTools
# print(C5_cdmTools)
