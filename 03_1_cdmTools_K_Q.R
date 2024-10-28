
# Estimate K and Q --------------------------------------------------------

# We have different packages available to estimate the K (the number of attribute)
# and Q-matrix 

# Introduction to packages --------------------------------------------------------

# {cdmTools}: developed in 2021 by Nájera, P.; Sorrel, M.A.; Abad, F.J https://cran.r-project.org/web/packages/cdmTools/cdmTools.pdf

# machine learning to estimate K 
# then using code developed in 2019 by Menta Chuang   
#  A Gibbs Sampling Algorithm that Estimates the Q-matrix for the DINA Model
# relate with paper in 2019 by Mengta Chuang https://www.sciencedirect.com/science/article/pii/S0022249619300264



# {cdmTools} ----------------------------------------------------

# Set a seed for reproducibility of random operations
set.seed(123)

cdmTools_K_Q_estimate <- function(cluster_subset) {
  # Ensure the data is in the correct matrix format
  data_matrix <- as.matrix(cluster_subset[,-1])
  
  # Estimate K using parallel analysis
  res_paK <- cdmTools::paK(data_matrix, cor = "cor")
  suggested_k <- res_paK$sug.K
  
  # Compare model fit across different values of K
  res_modelcompK <- cdmTools::modelcompK(data_matrix, exploreK = 1:6)
  
  # Plot BIC
  plot_bic <- ggplot(data = data.frame(K = 1:6, BIC = res_modelcompK$fit$BIC),
                     aes(x = K, y = BIC)) +
    geom_line() +
    geom_point() +
    labs(title = "BIC by Number of Clusters", x = "Number of Clusters", y = "BIC")
  
  # Plot AIC
  plot_aic <- ggplot(data = data.frame(K = 1:6, AIC = res_modelcompK$fit$AIC),
                     aes(x = K, y = AIC)) +
    geom_line() +
    geom_point() +
    labs(title = "AIC by Number of Clusters", x = "Number of Clusters", y = "AIC")
  
  # Extract and return the used Q matrix for the suggested K
  used_q <- res_modelcompK$usedQ[[paste("K", suggested_k, sep = "")]]
  
  # Visualizing the Q matrix
  plot_q <- as.data.frame(used_q)
  rownames(plot_q) <- c(paste0("item", 1:nrow(plot_q)))
  colnames(plot_q) <- c(paste0("attribute", 1:ncol(plot_q)))
  plot_q$Item <- rownames(plot_q)
  plot_q$Item <- factor(plot_q$Item, levels = (unique(plot_q$Item)))
  
  plot1 <- plot_q %>%
    tidyr::pivot_longer(cols = 1:ncol(used_q), names_to = "Attribute", values_to = "Accuracy")
  
  q_matrix_plot <- ggplot(plot1, aes(x = Attribute, y = Item)) +
    geom_tile(aes(fill = as.factor(Accuracy)), color = "white") +
    scale_fill_manual(values = c("0" = "white", "1" = "#0F68A8"),
                      labels = c("0" = "0", "1" = "1")) +
    labs(title = "Visualization of Q Matrix",
         fill = "Accuracy") +  # Set legend title
    theme_minimal() +
    labs(title = "Visualization of Q Matrix")
  
  # Create a descriptive summary
  output_description <- paste("The analysis suggests K =", suggested_k,
                              "as the optimal number of latent attributes based on parallel analysis.",
                              "Model comparisons for K=1 to K=6 are provided with respective AIC and BIC values plotted.",
                              "The used Q matrix for the suggested K is also returned.")
  
  # Return all relevant outputs in a list
  return(list(
    description = output_description,
    suggested_k = suggested_k,
    model_comparisons = res_modelcompK,
    plot_bic = plot_bic,
    plot_aic = plot_aic,
    used_q = used_q,
    q_matrix_plot = q_matrix_plot
  ))
}

# For example, cluster_subset is C5_subset
C5_cdmTools <- cdmTools_K_Q_estimate(C5_subset)

print(C5_cdmTools)



