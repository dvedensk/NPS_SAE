
functions {
  real bern_glm_ll_slice(array[] int y_slice,
                         int start, 
                         int end,
                         matrix X,
                         vector beta,
                         array[] int puma_id,
                          vector a_puma
                         ) {
    // intercept = 0 because X is built with ~ 0 + AGEP + SEX + RAC1P
    //return bernoulli_logit_glm_lpmf(y_slice | X[start:end, ], 0, beta); //bernoulli_logit_lpmf not much difference
  
  int N_slice = end - start + 1;
    vector[N_slice] eta = X[start:end, ] * beta;
    for (i in 1:N_slice) {
      eta[i] += a_puma[puma_id[start + i - 1]];
    }
    return bernoulli_logit_lpmf(y_slice | eta);
  }
  
  }


data {
  int<lower=1> n;
  int<lower=1> k;
  matrix[n, k] X;
  array[n] int<lower=0, upper=1> y;
  int<lower=1> grainsize;

  // prior scales for each coefficient
 // vector<lower=0>[k] beta_scale;
  
  
    // PUMA indexing
  int<lower=1> J;                 // number of PUMA levels
  array[n] int<lower=1, upper=J> puma_id;
  
  
  

}

parameters {
  vector[k] beta;
  // random intercept
  vector[J] z_puma;
  real<lower=0> sigma_puma;
}

transformed parameters {
  vector[J] a_puma = sigma_puma * z_puma;  // PUMA random effect
}

model {
  beta ~ normal(0, 3);
  z_puma ~ normal(0, 1);
 sigma_puma ~ cauchy(0, 5);  
 
  // Likelihood
  target += reduce_sum(bern_glm_ll_slice,
    y, grainsize,
    X, beta, puma_id, a_puma);
}

