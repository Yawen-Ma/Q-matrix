// =============================================================================
// A unified framework integrating Q matrix + DINA + latent transition model
// =============================================================================

// -----------------------------------------------------------------------------
// Attributes: 2 | Time points: 2 | Outcomes: 2 (binary: non-mastery, mastery)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// FUNCTIONS
// -----------------------------------------------------------------------------
functions {
  /**
   * Map 4 latent class 00, 10, 01, 11 (2 attributes with binary outcomes) to binary attribute vector of length 2
   */
  array[] int class2attr(int c) {
    array[2] int a;
    if (c == 1)        a = {0, 0}; // latent class 1: 00
    else if (c == 2)   a = {1, 0}; // latent class 2: 10
    else if (c == 3)   a = {0, 1}; // latent class 3: 01
    else               a = {1, 1}; // latent class 4: 11
    return a;
  }

  /**
   * calc_xi:
   * Compute the ideal response 0/1 for a given latent attribute vector (alpha)
   * and item Q-matrix row (Qrow) under a conjunctive (DINA) rule.
   * Returns:
   *     - If all required attributes are mastered (a[k]=1 for all required k),
   *       then ideal response = 1.
   *     - If any required attribute is not mastered (a[k]=0),
   *       then ideal response = 0.
   *     - If attribute is not required, it does not affect the response.
   */
  int calc_xi(array[] int a, array[] real Qrow) { 
    return ( (Qrow[1] > 0.5 ? a[1] : 1)
           * (Qrow[2] > 0.5 ? a[2] : 1) );
  }
}

// -----------------------------------------------------------------------------
// Data
// -----------------------------------------------------------------------------
data {
  int<lower=1> N;           // Number of examinees
  int<lower=1> J;           // Items per assessment
  int<lower=1> K;           // Number of attributes
  int<lower=1> P;           // Number of covariates

  array[N, J] int <lower=0,upper=1> Y1;  // Responses at time 1
  array[N, J] int <lower=0,upper=1> Y2;  // Responses at time 2
  matrix[N, P] Z;                // Covariate matrix
}

// -----------------------------------------------------------------------------
// Q matrix 
// + 
// Latent transition model: initial state parameters (beta) + transition parameters (gamma)
// -----------------------------------------------------------------------------
parameters {
  real<lower=0, upper=1> theta;                   // Hyperparameter for Q: sparsity
  array[J, K] real<lower=0, upper=1> Q_time1_raw;        // Raw Q‐matrix at time 1
  array[J, K] real<lower=0, upper=1> Q_time2_raw;        // Raw Q‐matrix at time 2

  vector<lower=0.05, upper=0.2>[J] g1;            // Guess parameters at time 1
  vector<lower=0.05, upper=0.2>[J] s1;            // Slip parameters at time 1
  vector<lower=0.05, upper=0.2>[J] g2;            // Guess parameters at time 2
  vector<lower=0.05, upper=0.2>[J] s2;            // Slip parameters at time 2

  vector[2] beta0;           // Intercepts for initial attribute probabilities
  matrix[2, P] betaZ;        // Coefficients for initial attribute probabilities

  matrix[2, P+1] gamma01;    // Transition prob (0->1) non-mastery to mastery
  matrix[2, P+1] gamma10;    // Transition prob (1->0) mastery to non-mastery
}

// -----------------------------------------------------------------------------
// TRANSFORMED PARAMETERS
// -----------------------------------------------------------------------------
transformed parameters {
  // 1) Constrain Q‐matrices
  matrix[J, K] Q_time1 = to_matrix(Q_time1_raw);
  matrix[J, K] Q_time2 = to_matrix(Q_time2_raw);
  // Time 1 fixed entries
  Q_time1[1,] = [1, 0];
  Q_time1[4,] = [0, 1];
  // Time 2 fixed entries
  Q_time2[1,] = [0, 1];
  Q_time2[2,] = [1, 0];

  // 2) Ideal responses eta at two time points
  // compute ideal response for each item j and latent class c at time 1 and 2
  array[J, 4] real xi1;
  array[J, 4] real xi2;
  for (j in 1:J){
    for (c in 1:4) {
      xi1[j, c] = calc_xi(class2attr(c), to_array_1d(Q_time1[j])); //xi1: for each item x for each latent class
      xi2[j, c] = calc_xi(class2attr(c), to_array_1d(Q_time2[j]));
    }
  }
  // 3) Attribute‐level probabilities
  matrix[N, 2] p1;     // Initial state probability
  matrix[N, 2] p01;    // Transiion probability (0->1)
  matrix[N, 2] p10;    // Transiion probability (1->0)
  for (n in 1:N) {
    for (k in 1:2) {
      // logistic regression: latent transition model
      p1[n, k]  = inv_logit(beta0[k] + dot_product(betaZ[k], Z[n])); 
      p01[n, k] = inv_logit(gamma01[k, 1] + dot_product(gamma01[k, 2:], Z[n]));
      p10[n, k] = inv_logit(gamma10[k, 1] + dot_product(gamma10[k, 2:], Z[n]));
    }
  }

  // 4) Latent class prior at time 1
  // log prior probability of student n at latent class c
  matrix[N, 4] log_nu1;
  for (n in 1:N) {
    real p11 = p1[n,1];
    real p12 = p1[n,2];
    log_nu1[n,1] = log1m(p11) + log1m(p12); // class 1: 00
    log_nu1[n,2] = log(p11)  + log1m(p12); // class 2: 10
    log_nu1[n,3] = log1m(p11) + log(p12); // class 3: 01
    log_nu1[n,4] = log(p11)  + log(p12); // class 4: 11 
  }

  // 5) 4×4 log transition matrix
  array[N, 4, 4] real log_trans;
  for (n in 1:N){
    for (c1 in 1:4) {
      array[2] int a1 = class2attr(c1); 
      for (c2 in 1:4) {
        array[2] int a2 = class2attr(c2); 
        real lp = 0;
        for (k in 1:2) {
          if (a1[k]==0 && a2[k]==0) lp += log1m(p01[n,k]);       // 0->0: x(1-p01)
          else if (a1[k]==0 && a2[k]==1) lp += log(p01[n,k]);    // 0->1: xp01
          else if (a1[k]==1 && a2[k]==1) lp += log1m(p10[n,k]);  // 1->1: x(1-p10)
          else                           lp += log(p10[n,k]);    // 1->0: xp10
        }
        log_trans[n,c1,c2] = lp;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// MODEL
// -----------------------------------------------------------------------------
model {
  // Priors
  theta ~ beta(3, 1);
  for (j in 1:J) {
    if (j != 1 && j != 4)
      Q_time1_raw[j] ~ beta(theta + 1e-3, (1-theta) + 1e-3);
    if (j != 1 && j != 2)
      Q_time2_raw[j] ~ beta(theta + 1e-3, (1-theta) + 1e-3);
  }
  beta0   ~ normal(0, 0.6);
  to_vector(betaZ)  ~ normal(0, 0.3);
  gamma01[,1] ~ normal(-1.5, 0.3);
  gamma10[,1] ~ normal(-3.0, 1.0);
  to_vector(gamma01[,2:]) ~ normal(0, 0.3);
  to_vector(gamma10[,2:]) ~ normal(0, 0.3);

  // Likelihood
  for (n in 1:N) {
    vector[4] ll1;
    for (c1 in 1:4) {
      real lp1 = log_nu1[n,c1];
      for (j in 1:J){
        lp1 += bernoulli_lpmf(Y1[n,j] | xi1[j,c1]*(1-s1[j]) + (1-xi1[j,c1])*g1[j]); 
        }
      vector[4] ll2;
      for (c2 in 1:4) {
        real lp2 = log_trans[n,c1,c2];
        for (j in 1:J){
          lp2 += bernoulli_lpmf(Y2[n,j] | xi2[j,c2]*(1-s2[j]) + (1-xi2[j,c2])*g2[j]);
        }
        ll2[c2] = lp2;
      }
      ll1[c1] = lp1 + log_sum_exp(ll2);
    }
    target += log_sum_exp(ll1);
  }
}

// -----------------------------------------------------------------------------
// GENERATED QUANTITIES
// -----------------------------------------------------------------------------
generated quantities {
  matrix[N, 4] prob_class1;
  matrix[N, 4] prob_class2;
  matrix[N, 2] prob_attr1;
  matrix[N, 2] prob_attr2;

  for (n in 1:N) {
    // Recompute log‐likelihoods for posterior
    vector[4] lp1;
    array[4, 4] real lp2_arr; 
    for (c1 in 1:4) {
      lp1[c1] = log_nu1[n,c1];
      for (j in 1:J) {
        lp1[c1] += bernoulli_lpmf(Y1[n,j] | xi1[j,c1]*(1-s1[j]) + (1-xi1[j,c1])*g1[j]);
      }
      for (c2 in 1:4) {
        lp2_arr[c1,c2] = log_trans[n,c1,c2];
        for (j in 1:J){
          lp2_arr[c1,c2] += bernoulli_lpmf(Y2[n,j] | xi2[j,c2]*(1-s2[j]) + (1-xi2[j,c2])*g2[j]);
        }
      }
    }

    // Posterior over latent classes at time 1
    prob_class1[n] = softmax(lp1)';

    // Posterior over latent classes at time 2
    vector[4] num2 = rep_vector(-1e30, 4);
    for (c1 in 1:4) {
      for (c2 in 1:4) {
        num2[c2] = log_sum_exp(num2[c2], lp2_arr[c1,c2] + lp1[c1]);
        }
    }
    prob_class2[n] = softmax(num2)';

    // Marginal attribute posteriors
    for (k in 1:2) {
      real a1 = 0; real a2 = 0;
      for (c in 1:4) {
        a1 += prob_class1[n,c] * class2attr(c)[k];
        a2 += prob_class2[n,c] * class2attr(c)[k];
      }
      prob_attr1[n,k] = a1;
      prob_attr2[n,k] = a2;
    }
  }
}






