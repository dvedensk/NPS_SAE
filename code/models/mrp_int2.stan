
functions {
  real bern_ll_slice(array[] int y_slice, int start, int end,
                     matrix X,
                     array[] int psi_bin,
                     vector lp_psi,
                     vector beta,
                     real beta_psi,
                     vector zeta,
                     array[] int puma_id,
                     vector a_puma) {
    int M = end - start + 1;
    vector[M] alpha_slice =
        beta_psi * lp_psi[start:end]
      + zeta[psi_bin[start:end]]
      + a_puma[puma_id[start:end]];

    return bernoulli_logit_glm_lpmf(y_slice | X[start:end, ], alpha_slice, beta);
 
  }
}

data {
  // --- Training data 
  int<lower=1> n;
  int<lower=1> k;
  matrix[n, k] X;
  array[n] int<lower=0, upper=1> y;
  int<lower=1> J;
  array[n] int<lower=1, upper=J> puma_id;
  int<lower=1> G;    
  vector[n] lp_psi; 
  // number of psi bins
  array[n] int<lower=1, upper=G> psi_bin;   // bin index for each training unit
  // reduce_sum grainsize
  int<lower=1> grainsize;
  // ---prior
 // vector[k] beta_scale;      // Vector of scales for beta
//  real beta_psi_scale;       // Scale for beta_psi
  real sigma_psi_rate;       // Rate for sigma_psi
  
}


parameters {
  vector[k] beta;
  real beta_psi;                
//real<lower=0, upper=5> sigma_psi;
  vector[G]     zeta_raw;
  real<lower=0> sigma_psi;
    real<lower=0> sigma_puma;

vector[J] z_puma;

}

transformed parameters {
  vector[G] zeta = sigma_psi * zeta_raw;    
  vector[J] a_puma = sigma_puma * z_puma;
}

model {
  beta      ~ normal(0, 3);
  beta_psi  ~ normal(0, 3);
  sigma_psi ~ cauchy(0, 5);  
  zeta_raw  ~ normal(0, 1);
  z_puma     ~ normal(0, 1);
 sigma_puma ~ cauchy(0, 5);  
  target += reduce_sum(bern_ll_slice, y, grainsize,
                       X, psi_bin, lp_psi, beta, beta_psi, zeta,
                       puma_id, a_puma);
                       
}



