functions {
  real bern_ll_slice(array[] int y_slice, int start, int end,
                     matrix X,
                     array[] int psi_bin,
                     vector lp_psi,
                     vector beta,
                     real beta_psi,
                     vector zeta) {
    int M = end - start + 1;
    // vector[M] eta;
    // // fixed effects + selection term + psi-bin 
    // eta =  X[start:end, ] * beta
    //      + beta_psi * lp_psi[start:end]
    //      + zeta[psi_bin[start:end]];
    // return bernoulli_logit_lpmf(y_slice | eta);
    
    vector[M] alpha_slice = beta_psi * lp_psi[start:end] + zeta[psi_bin[start:end]];
    return bernoulli_logit_glm_lpmf(y_slice | X[start:end, ], alpha_slice, beta);
  }
}
data {
  // --- Training data 
  int<lower=1> n;
  int<lower=1> k;
  matrix[n, k] X;
  array[n] int<lower=0, upper=1> y;

  // inclusion prob & bin for MRP-INT
  // vector<lower=0, upper=1>[n] psi;          // estimated inclusion probabilities for training units
  int<lower=1> G;    
  vector[n] lp_psi; 
  // number of psi bins
  array[n] int<lower=1, upper=G> psi_bin;   // bin index for each training unit

  // reduce_sum grainsize
  int<lower=1> grainsize;
  // ---prior
  vector[k] beta_scale;      // Vector of scales for beta
  real beta_psi_scale;       // Scale for beta_psi
  real sigma_psi_rate;       // Rate for sigma_psi
}


parameters {
  vector[k] beta;
  real beta_psi;                    // fixed effect for logit(psi)
  // psi-bin random effects
real<lower=0, upper=5> sigma_psi;
  vector[G]     zeta_raw;
}

transformed parameters {
  vector[G] zeta = sigma_psi * zeta_raw;    // psi-bin RE
}


model {
  // Data-informed priors
  beta      ~ normal(0, beta_scale);         
  beta_psi  ~ normal(0, beta_psi_scale);
  sigma_psi ~ exponential(sigma_psi_rate);    
  zeta_raw  ~ normal(0, 1);

  // Parallelized likelihood
  target += reduce_sum(bern_ll_slice, y, grainsize,
                       X, psi_bin, lp_psi, beta, beta_psi, zeta);
}