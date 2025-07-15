// =============================================================================
//  DINA, 2 Attribute, 4 Class, 2 Time Points, Covariates Implementation
// =============================================================================
//  – K = 2  gives  C = 4 latent classes (00,10,01,11)
//  – Discrete attribute patterns; likelihood marginalised analytically.
//  – Initial mastery (beta) and transitions (gamma) depend on covariates Z.
//  – Slip/guess item parameters per occasion with Uniform(.05,.15) priors.
//  – Generated quantities return posterior probabilities for classes & skills.
// -----------------------------------------------------------------------------
//  NOTE:  Works **only for K = 2**.  For more attributes use the generic file.
// =============================================================================
functions {
  array[] int class2attr(int c) {
    array[2] int a;
    if      (c==1) a={0,0};
    else if (c==2) a={1,0};
    else if (c==3) a={0,1};
    else           a={1,1};
    return a;
  }
  int calc_xi(array[] int a, array[] real Qrow){
    return ( (Qrow[1]>0.5 ? a[1]:1) * (Qrow[2]>0.5 ? a[2]:1) ); // if true, take [1], otherwise 1
  }
}


// =============================================================================
// DATA BLOCK
// =============================================================================

data {
  int<lower=1> N;                       // students
  int<lower=1> J;                       // items per occasion
  int<lower=1> P;                       // covariates

  array[N,J] int<lower=0,upper=1> Y1;   // responses time‑1
  array[N,J] int<lower=0,upper=1> Y2;   // responses time‑2
  matrix[N,P] Z;                        // covariate matrix
  
  array[J,2] real<lower=0,upper=1> Q1;   // Q-matrix time-1
  array[J,2] real<lower=0,upper=1> Q2;   // Q-matrix time-2
}

transformed data {
  /* ---------- 4. xi(i,c) for two times ---------- */
  array[J,4] int xi1;
  array[J,4] int xi2;

  for (i in 1:J) {
    for (c in 1:4) {
      xi1[i,c] = calc_xi(class2attr(c), Q1[i]);
      xi2[i,c] = calc_xi(class2attr(c), Q2[i]);
    }
  }
}

// =============================================================================
// PARAMETERS
// =============================================================================

parameters {
 array[J]  real<lower = 0.05, upper = 0.15> g1;
 array[J]  real<lower = 0.05, upper = 0.15> s1;
 array[J]  real<lower = 0.05, upper = 0.15> g2;
 array[J]  real<lower = 0.05, upper = 0.15> s2;
// regressions
  vector[2] beta0;
  matrix[2,P] betaZ;
  matrix[2,P+1] gamma01;
  matrix[2,P+1] gamma10;
}

// =============================================================================
// TRANSFORMED PARAMETERS
// =============================================================================

transformed parameters {
  // ---- Attribute‑level probabilities --------------------------------------
  matrix[N,2] p1;
  matrix[N,2] p01;
  matrix[N,2] p10;
  
  for (n in 1:N) {
    for (k in 1:2) {
      p1 [n,k] = inv_logit( beta0[k] + dot_product(betaZ[k] , Z[n]) );
      p01[n,k] = inv_logit( gamma01[k,1] + dot_product(gamma01[k,2:(P+1)], Z[n]) );
      p10[n,k] = inv_logit( gamma10[k,1] + dot_product(gamma10[k,2:(P+1)], Z[n]) );
    }
  }

  // ---- Time‑1 class log prior ---------------------------------------------
  matrix[N,4] log_nu1;
  for (n in 1:N) {
    real p11 = p1[n,1];
    real p12 = p1[n,2];
    log_nu1[n,1] = log1m(p11) + log1m(p12);
    log_nu1[n,2] = log(p11)  + log1m(p12);
    log_nu1[n,3] = log1m(p11)+ log(p12);
    log_nu1[n,4] = log(p11)  + log(p12);
  }

  // ---- 4×4 log transition --------------------------------------------------
  array[N,4,4] real log_trans;
  for (n in 1:N) {
    for (c1 in 1:4) {
      array[2] int a1 = class2attr(c1);
      for (c2 in 1:4) {
        array[2] int a2 = class2attr(c2);
        real lp = 0;
        for (k in 1:2) {
          if      (a1[k]==0 && a2[k]==0) lp += log1m(p01[n,k]);
          else if (a1[k]==0 && a2[k]==1) lp += log(p01[n,k]);
          else if (a1[k]==1 && a2[k]==1) lp += log1m(p10[n,k]);
          else                           lp += log(p10[n,k]);
        }
        log_trans[n,c1,c2] = lp;
      }
    }
  }
}

// =============================================================================
// MODEL
// =============================================================================

model {
  // ---- Priors --------------------------------------------------------------
  beta0            ~ normal(0,0.6);
  to_vector(betaZ) ~ normal(0,0.3);
  gamma01[,1] ~ normal(-1.5, 0.3);   // learning
  gamma10[,1] ~ normal(-3, 1.0);   // forgetting
  to_vector(gamma01[, 2:(P+1)]) ~ normal(0, 0.3);
  to_vector(gamma10[, 2:(P+1)]) ~ normal(0, 0.3);
 
       
  // ---- Likelihood ----------------------------------------------------------
  for (n in 1:N) {
    vector[4] log_like1;

    for (c1 in 1:4) {
      real lp1 = log_nu1[n,c1];
      for (i in 1:J) {
        real p = xi1[i,c1]*(1-s1[i]) + (1-xi1[i,c1])*g1[i];
        lp1 += bernoulli_lpmf(Y1[n,i] | p);
      }

      vector[4] log_like2;
      for (c2 in 1:4) {
        real lp2 = log_trans[n,c1,c2];
        for (i in 1:J) {
          real p = xi2[i,c2]*(1-s2[i]) + (1-xi2[i,c2])*g2[i];
          lp2     += bernoulli_lpmf(Y2[n,i] | p);
        }
        log_like2[c2] = lp2;
      }
       log_like1[c1] = lp1 + log_sum_exp(log_like2);
    }
    target += log_sum_exp(log_like1);
    }
  }

// =============================================================================
// GENERATED QUANTITIES
// =============================================================================

generated quantities {
  matrix[N,4] prob_class1;
  matrix[N,4] prob_class2;
  matrix[N,2] prob_attr1;
  matrix[N,2] prob_attr2;

  for (n in 1:N) {
    vector[4] log_like1;
    array[4] vector[4] log_like2;   
    vector[4] lp1_arr;              

    for (c1 in 1:4) {
      real lp1 = log_nu1[n,c1];
      for (i in 1:J) {
        real p = xi1[i,c1]*(1-s1[i]) + (1-xi1[i,c1])*g1[i];
        lp1 += bernoulli_lpmf(Y1[n,i] | p);
      }
      lp1_arr[c1] = lp1;            // ★ 存起来

      for (c2 in 1:4) {
        real lp2 = log_trans[n,c1,c2];
        for (i in 1:J) {
          real p = xi2[i,c2]*(1-s2[i]) + (1-xi2[i,c2])*g2[i];
          lp2 += bernoulli_lpmf(Y2[n,i] | p);
        }
        log_like2[c1][c2] = lp2;
      }
      log_like1[c1] = lp1 + log_sum_exp(log_like2[c1]);
    }

    /* ---- L1 posterior ---- */
    prob_class1[n] = softmax(log_like1)';
    /* ---- L2 posterior ---- */
    vector[4] num_l2 = rep_vector(-1e30, 4);  // log numerators
    for (c1 in 1:4) {
      for (c2 in 1:4) {
      num_l2[c2] = log_sum_exp(
                    num_l2[c2],
                    log_like2[c1][c2] + lp1_arr[c1]   
                 );
          }
    }
    prob_class2[n] = softmax(num_l2)';

    /* ---- draw prob_attr1 and prob_attr2 as posterior ---- */
    for (k in 1:2) {
      real a1 = 0;
      real a2 = 0;
      for (c in 1:4) {
        a1 += prob_class1[n,c] * class2attr(c)[k];
        a2 += prob_class2[n,c] * class2attr(c)[k];
      }
      prob_attr1[n,k] = a1;
      prob_attr2[n,k] = a2;
    }
  }
}


