data {
  int<lower=1> J; //num cells
  array[J] int<lower=0> n; //sample sizes by cell
  vector<lower=0>[J] N; // population cell sizes
  int<lower=0> K;  //num covariates
  matrix[J, K] X; 

  int<lower=1> G_puma; //number of pumas
  array[J] int<lower=1, upper=G_puma> puma_id;
}

parameters {
  real alpha0;
  vector[K] alpha;

  vector[G_puma] z_puma;
  real<lower=0> sigma_puma;
}

transformed parameters {
  vector[G_puma] b_puma;
  vector[J] eta;
  vector<lower=0, upper=1>[J] psi; //propensity for each cell
  simplex[J] p;

  b_puma = sigma_puma * z_puma;
  eta = alpha0 + X * alpha + b_puma[puma_id];

  psi = inv_logit(eta);
  //p = (N .* psi) / dot_product(N, psi);
  p = softmax(log(N) + log_inv_logit(eta)); //same as above but more stable

}

model {
  alpha0 ~ normal(0, 3);
  alpha ~ normal(0, 3);

  z_puma ~ normal(0,1);
  sigma_puma ~ cauchy(0, 5);

  n ~ multinomial(p);
}

generated quantities {
  vector[J] logit_psi;
  logit_psi = eta;
}

