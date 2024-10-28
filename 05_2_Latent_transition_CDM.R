
# Latent Transition CDMs --------------------------------------------------

## Define Time point ------------------------------------------------------

yearly_points <- define_time_point(items_C5, "2021-08-01", "2023-08-01", 12)

## Analysis ------------------------------------------------------

library(GDINA)


Year1 <- yearly_points %>% 
  filter(Time_Period == "Year 1") %>%
  group_by(USER_ID, skill_family) %>%
  mutate(
         response_accuracy = as.integer(is_level_mastered)
        ) %>%
  ungroup() %>%
  # Take out response accuracy at this stage for generating q-matrix
  select(USER_ID, Item_at_time_point, response_accuracy) %>%
  # Restructure the dataset, user by items (32 columns for each item)
  pivot_wider(names_from = Item_at_time_point, values_from = response_accuracy, id_cols = USER_ID) %>%
  # Pick up the first 9 responses from year 1
  .[, 1:10] %>%
  na.omit()

Q_estimate <- cdmTools::modelcompK(Year1[, -1], exploreK = 1:4)

# Fit the response data at Time point 1 to the selected models
fit.t1 <- GDINA(dat = Year1[, -1], Q = Q_estimate$usedQ$K4, model = "GDINA", mono.constraint = TRUE, verbose = 0)


Year2 <- yearly_points %>% 
  filter(USER_ID %in% Year1$USER_ID) %>%
  filter(Time_Period == "Year 2") %>%
  group_by(USER_ID, skill_family) %>%
  mutate(
    response_accuracy = as.integer(is_level_mastered)
  ) %>%
  ungroup() %>%
  # Take out response accuracy at this stage for generating q-matrix
  select(USER_ID, Item_at_time_point, response_accuracy) %>%
  # Restructure the dataset, user by items (32 columns for each item)
  pivot_wider(names_from = Item_at_time_point, values_from = response_accuracy, id_cols = USER_ID) %>%
  # Pick up the first 9 responses from year 2
  .[, 1:10] %>%
  na.omit()

Q_estimate_T2 <- cdmTools::modelcompK(Year2[, -1], exploreK = 1:4)

# Fit the response data at Time point 2 to the selected models
fit.t2 <- GDINA(dat = Year2[, -1], Q = Q_estimate_T2$usedQ$K4, model = "GDINA", mono.constraint = TRUE, verbose = 0)

fit.object = list()
fit.object[[1]] <- fit.t1
fit.object[[2]] <- fit.t2

t = 2 # the number of time points
K = 4 # the number of attributes
N = nrow(Year2) # the number of observations
cep = CEP_t(fit.object = fit.object, t = t, K = K, N = N)

# The CEP matrices of the attributes
cep$cep.matrix


####################################################################
## ::::::::::::::::: From package, not work ::::::::::: ###
####################################################################


z_t1_test = matrix(personparm(fit.t1)[1:1000,], ncol = 4)
z_t2_test = matrix(personparm(fit.t2)[1:1000,], ncol = 4)
#Initial values of initial state's regression coefficients
beta_in = matrix(0, ncol(z_t1_test), K)

# Initial values of transition probability's regression coefficients
# These were computed using the raw data.
# When Gender coding is 1 = male, 0 = female:
ga01_in = cbind(c(-2.15, 0.56, 0.09, -0.79),
                c(-1.6, 0.05, -0.01, -0.38),
                c(-1.25, 0.06, -0.25, 0.14),
                c(-1.18, -0.26, 0.04, 0.37),
                )
#initial values of regression coefficients (for transition from 0 to 1)
ga10_in = cbind(c(-0.84, -0.18, -0.14, 0.23),
                c(-0.18, 0.49, 0.44, -0.35),
                c(-0.22, 0.18, 0.37, -0.45),
                c(-0.49, 0.10, 0.43, 0.20))  

step3.output_test <- step3.est(cep = cep, z_t1 = z_t1_test, z_t2 = z_t2_test,
                               K = K, t = t, beta_in, ga01_in, ga10_in)


###############################################################################
## :::::::::::::::::  {LMest} to get the transition probability ::::::::::: ###
###############################################################################

# take 1000 users randomly
z_t1_test = matrix(personparm(fit.t1)[1:1000,], ncol = 4)
z_t2_test = matrix(personparm(fit.t2)[1:1000,], ncol = 4)

library(LMest)

# Assuming data is already loaded and structured correctly
# Example: data = data.frame(student_id = c(1,1,2,2), time = c(1,2,1,2), skill1 = c(1,1,0,1), skill2 = c(1,0,1,1), skill3 = c(1,1,0,0), skill4 = c(1,0,1,1))
long_data$
fit$transition_probs

