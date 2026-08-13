data {
  int<lower=1> J;
  int<lower=1> J_obs;
  array[J_obs] int<lower=1, upper=J> obs_id;

  array[J_obs] int<lower=0> y;
  array[J_obs] int<lower=0> m;

  int<lower=0> K;
  matrix[J, K] X;
  vector[J] logit_psi;

  int<lower=1> G_puma;  //number of PUMAs
  array[J] int<lower=1, upper=G_puma> puma_id;
 
  int<lower=1> G_psi; // number of propensity score bins
  array[J] int<lower=1, upper=G_psi> psi_group_id;

  vector<lower=0>[J] N; //population cell sizes

  matrix[G_puma, J] A_puma; //matrix to map cells to PUMAs for post stratification 
}

parameters {
  real alpha0;
  vector[K] beta;
  real gamma;

  vector[G_puma] z_puma;
  real<lower=0> sigma_puma;

  vector[G_psi] z_psi_group;
  real<lower=0> sigma_psi_group;
}

transformed parameters {
  vector[G_puma] b_puma; //puma level random effects
  vector[G_psi] b_psi_group; // binned propensity score random effects
  vector[J] eta; // linear predictor

  b_puma = sigma_puma * z_puma;
  b_psi_group = sigma_psi_group * z_psi_group;

  eta = alpha0 +
    X * beta +
    gamma * logit_psi +
    b_puma[puma_id] +
    b_psi_group[psi_group_id];
}

model {
  alpha0 ~ normal(0, 3);
  beta ~ normal(0, 3);
  gamma ~ normal(0, 3);

  sigma_puma ~ cauchy(0, 5);
  z_puma ~ normal(0, 1);

  sigma_psi_group ~ cauchy(0, 5);
  z_psi_group ~ normal(0, 1);

  y ~ binomial_logit(m, eta[obs_id]);
}

generated quantities {
  //generate poststratified estimates at the PUMA level
  vector[G_puma] ybar_puma; 

  {
    vector[J] theta;
    vector[G_puma] numer;
    vector[G_puma] denom;

    theta = inv_logit(eta);

    numer = A_puma * (N .* theta);
    denom = A_puma * N;

    ybar_puma = numer ./ denom;
  }
}
