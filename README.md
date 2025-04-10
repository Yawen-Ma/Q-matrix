# RUN the most recent two files named "50rep_200N_6J_0.5theta.R" and "20_stan_0.5.stan".



Q-matrix project
│   Q-matrix.Rproj
│   README.md
└───data
│   │   C1.csv
│   │   C2.csv
│   │   C3.csv
│   │   C4.csv
│   │   C5.csv
│   │   C6.csv
│   │   C7.csv
│   │   C8.csv
│   │   C9.csv
└───plots
│   │   residuals.png
│   │   outcome_by_age.png
│   │   outcome_by_occupation.png
└───R
│   │   00_packages.R
│   │   01_load_data.R
│   │   02_items_response.R
│   │   03.1_cdmTools_K_Q.R
│   │   03.2_edina_K_Q.R
│   │   03.3_ML_.R (developing)
│   │   04_validate.R (developing)
│   │   05_CDM_model.R
│   │   06_final_plots_tables.R


# Estimate K and Q --------------------------------------------------------

# We have different packages available to estimate the K (the number of attribute)
# and Q-matrix 

# Introduction to packages --------------------------------------------------------

# {edina}: developed in 2020 by James Joseph Balamuta, Steven Andrew Culpepper, Jeffrey A. Douglas https://cran.r-project.org/web/packages/edina/edina.pdf
# relate with paper in 2018 "Bayesian estimation of the DINA Q-matrix" (Y. Chen, S.A. Culpepper, Y. Chen, J.A. Douglas) see https://link.springer.com/article/10.1007/s11336-017-9579-4

# {cdmTools}: developed in 2021 by Nájera, P.; Sorrel, M.A.; Abad, F.J https://cran.r-project.org/web/packages/cdmTools/cdmTools.pdf

# machine learning to estimate K 
# then using code developed in 2019 by Menta Chuang   
# A Gibbs Sampling Algorithm that Estimates the Q-matrix for the DINA Model
# relate with paper in 2019 by Mengta Chuang https://www.sciencedirect.com/science/article/pii/S0022249619300264

