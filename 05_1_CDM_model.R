
# CDM for RT RA -----------------------------------------------------------

# The Q-matrix obtained earlier
Q_matrix = C5_cdmTools$used_q

# The number of users
n_examinees = n_distinct(C5_subset$USER_ID)
n_items = J

# Design array (N-by-J-by-L) 	
# A N-by-J-by-L array indicating whether item j is administered to examiee i at l time point.
Response <- as.matrix(C5_subset[, -1])  # Convert to matrix first
response <- array(Response, dim = c(dim(Response), 1)) # Convert matrix to array

uni.testorder <- matrix(c(1,1), ncol=2, byrow=T)
uni.testver <- c(rep(1, n_examinees))

C5_sub_DINA_FOHM <- hmcdm(
  response,
  Q_matrix = C5_cdmTools$used_q,
  Test_order = uni.testorder,
  Test_versions = uni.testver,
  model = "DINA_FOHM",
  chain_length = 300, # Number of MCMC iterations
  burn_in = 100 # Burn-in period
  )



