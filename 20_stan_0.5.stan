//

functions {
  real compute_eta_AK_K_soft(int K, array[] real alpha_vec, array[] real Q_row) {
    real log_eta = 0; 
    for(k in 1:K) {
      log_eta += Q_row[k] * log(alpha_vec[k] + 1e-12);
    }
    return exp(log_eta);
  }
}

data {
  int<lower=1> N;                 // number of students
  int<lower=1> J;                 // number of items
  int<lower=1> K;                 // number of attributes
  
  // NEW SYNTAX: array[N, J] int<...>
  array[N, J] int<lower=0, upper=1> Y1;  // responses at time 1
  array[N, J] int<lower=0, upper=1> Y2;  // responses at time 2
  
  int<lower=1> P;                 // number of covariates
  matrix[N, P] Z;                 // N x P covariates
}

parameters {
  real<lower=0, upper=1> theta;

  array[J, K] real<lower=0, upper=1> Q_time1_raw;
  array[J, K] real<lower=0, upper=1> Q_time2_raw;

  array[J] real<lower=0, upper=1> g1;
  array[J] real<lower=0, upper=1> s1;
  array[J] real<lower=0, upper=1> g2;
  array[J] real<lower=0, upper=1> s2;

   // initial state regression coefficients
  // alpha1[i,k] = logistic( beta0[k] + sum_p betaZ[k,p]*Z[i,p] )
  vector[K] beta0;        
  matrix[K, P] betaZ;     
  
  // transition regression coefficients
  // p_01[k, i] = logistic( gamma01[k,1] + sum_p gamma01[k,p+1]*Z[i,p] )
  // p_10[k, i] = logistic( gamma10[k,1] + sum_p gamma10[k,p+1]*Z[i,p] )
  matrix[K, P+1] gamma01;
  matrix[K, P+1] gamma10;
}

transformed parameters {
  
   // The final Q_time1 and Q_time2 used in the measurement model
  array[J, K] real<lower=0, upper=1> Q_time1;
  array[J, K] real<lower=0, upper=1> Q_time2;

  // --------------------------------------------------------------------
  // 1) alpha1[i,k]
  //    We define alpha1 via logistic regression on covariates Z.
  // --------------------------------------------------------------------
  array[N, K] real<lower=0, upper=1> alpha1;
  
  // 2) p_01[i,k], p_10[i,k]: "transition" for each attribute
  //    p_01[i,k] = Probability alpha2=1 given alpha1=0
  //    p_10[i,k] = Probability alpha2=1 given alpha1=1
  array[N, K] real<lower=0, upper=1> p_01;
  array[N, K] real<lower=0, upper=1> p_10;
  
  // 3) alpha2[i,k]
  //    alpha2[i,k] = alpha1[i,k]*(1 - p_10[i,k]) + (1-alpha1[i,k])*p_01[i,k]
  array[N, K] real<lower=0, upper=1> alpha2;
  
  // compute alpha1
  for(i in 1:N) {
    for(k_ in 1:K) {
      real lin = beta0[k_];
      for(p_ in 1:P) {
        lin += betaZ[k_, p_] * Z[i, p_];
      }
      alpha1[i, k_] = inv_logit(lin); 
    }
  }
  
  // compute p_01, p_10
  for(i in 1:N) {
    for(k_ in 1:K) {
      {
        real lin01 = gamma01[k_,1];
        real lin10 = gamma10[k_,1];
        for(pp in 1:P) {
          lin01 += gamma01[k_, pp+1] * Z[i, pp];
          lin10 += gamma10[k_, pp+1] * Z[i, pp];
        }
        p_01[i,k_] = inv_logit(lin01);
        p_10[i,k_] = inv_logit(lin10);
      }
    }
  }
  
  // compute alpha2
  for(i in 1:N) {
    for(k_ in 1:K) {
      // if alpha1=0 => alpha2=1 with prob p_01
      // if alpha1=1 => alpha2=0 with prob p_10
      alpha2[i, k_] = alpha1[i,k_]*(1 - p_10[i,k_]) + (1 - alpha1[i,k_])* p_01[i,k_];
    }
  }
  
  // 
  for (j in 1:J) {
    for (k_ in 1:K) {
      Q_time1[j, k_] = Q_time1_raw[j, k_];
      Q_time2[j, k_] = Q_time2_raw[j, k_];
    }
  }

  // Time 1 constraints:
  Q_time1[1,1] = 1;  // Q_TIME1[1,1] <- 1
  Q_time1[1,2] = 0;  // Q_TIME1[1,2] <- 0
  Q_time1[4,1] = 0;  // Q_TIME1[4,1] <- 0
  Q_time1[4,2] = 1;  // Q_TIME1[4,2] <- 1

  // Time 2 constraints:
  Q_time2[1,1] = 0;  // Q_TIME2[1,1] <- 0
  Q_time2[1,2] = 1;  // Q_TIME2[1,2] <- 1
  Q_time2[2,1] = 1;  // Q_TIME2[2,1] <- 1
  Q_time2[2,2] = 0;  // Q_TIME2[2,2] <- 0

}

model {
  // 1) Prior on theta
 theta ~ beta(2,2);
  
// 2) approximate Bernoulli => Beta( c*theta, c*(1-theta) )
  for(j in 1:J) {
    //for(k_ in 1:K) {
      // skip forced constraints for Time1 row1 => (1,0), row4 => (0,1)
      // skip forced constraints for Time2 row1 => (0,1), row2 => (1,0)

      // For Time1:
      if(!(
           (j == 1 ) ||
           (j == 4 )
         )) {
       Q_time1_raw[j] ~ beta(1*theta + 1e-3, 1*(1-theta) + 1e-3);
      }

      // For Time2:
      if(!(
           (j == 1 ) ||
           (j == 2 )
         )) {
       Q_time2_raw[j] ~ beta(1*theta + 1e-3, 1*(1-theta) + 1e-3);
       
    //  }
    }
  }

//-------------------------------------------------
  // 3) Other Priors (slip, guess, alpha, etc.)
  //-------------------------------------------------
  
  // A) item parameters
  for(j_ in 1:J) {
    g1[j_] ~ beta(5, 15);
    s1[j_] ~ beta(5, 15);
    g2[j_] ~ beta(5, 15);
    s2[j_] ~ beta(5, 15);
  }
  
  // B) logistic regression for alpha1
  for(k_ in 1:K) {
    beta0[k_] ~ normal(0, 0.3);
    for(p_ in 1:P) {
      betaZ[k_, p_] ~ normal(0, 0.3);
    }
  }
  
  // C) transition coefficients for p_01, p_10
  for(k_ in 1:K) {
    for(pp in 1:(P+1)) {
      gamma01[k_, pp] ~ normal(0, 0.3);
      gamma10[k_, pp] ~ normal(0, 0.3);
    }
  }
  
  
  //=================
  // B) Likelihood
  //=================
  // Y1 measurement model:
  //   prob of correct = eta1 * (1 - s1[j]) + (1-eta1)*g1[j]
  //   where eta1 = product_{k} alpha1[i,k]^Q_time1[j,k]
  
  for(i in 1:N) {
    for(j_ in 1:J) {
      // Soft Q + Soft alpha => use compute_eta_AK_K_soft
      real eta1 = compute_eta_AK_K_soft(K,
        to_array_1d(alpha1[i]),  // alpha1[i,*] in [0,1]
        to_array_1d(Q_time1[j_]) // Q_time1[j_,*] in [0,1]
      );
      real p_j = eta1*(1 - s1[j_]) + (1-eta1)*g1[j_];
      Y1[i, j_] ~ bernoulli(p_j);
    }
  }

  
  for(i in 1:N) {
    for(j_ in 1:J) {
      real eta2 = compute_eta_AK_K_soft(K,
        to_array_1d(alpha2[i]),
        to_array_1d(Q_time2[j_])
      );
      real p_j = eta2*(1 - s2[j_]) + (1-eta2)*g2[j_];
      Y2[i, j_] ~ bernoulli(p_j);
    }
  }
}


