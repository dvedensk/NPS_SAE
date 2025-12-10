
functions {
  real bern_glm_ll_slice(array[] int y_slice, int start, int end,
                         matrix X,
                         vector beta) {
    // intercept = 0 because X is built with ~ 0 + AGEP + SEX + RAC1P
    return bernoulli_logit_glm_lpmf(y_slice | X[start:end, ], 0, beta); //bernoulli_logit_lpmf not much difference
  }
}

data {
  int<lower=1> n;
  int<lower=1> k;
  matrix[n, k] X;
  array[n] int<lower=0, upper=1> y;
  int<lower=1> grainsize;

  // NEW: prior scales for each coefficient
  vector<lower=0>[k] beta_scale;
}

parameters {
  vector[k] beta;
}

model {
  // Priors: beta ~ Normal(0, beta_scale[j])
  beta ~ normal(0, beta_scale);

  // Likelihood
  target += reduce_sum(bern_glm_ll_slice, y, grainsize, X, beta);
}

