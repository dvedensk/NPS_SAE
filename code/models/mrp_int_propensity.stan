data {
  int<lower=1> J; //num cells
  array[J] int<lower=0> n; //sample sizes by cell
  vector<lower=0>[J] N; // population cell sizes
  int<lower=0> K;  //num covariates
  matrix[J, K] X; 
}

parameters {
  real alpha0; //intercept
  vector[K] alpha; 
}

transformed parameters {
  vector[J] eta;
  vector<lower=0, upper=1>[J] psi; //propensity for each cell
  simplex[J] p;

  eta = alpha0 + X * alpha;
  psi = inv_logit(eta);
  //p = (N .* psi) / dot_product(N, psi);
  p = softmax(log(N) + log_inv_logit(eta)); //same as above but more stable

}

model {
  alpha0 ~ normal(0, 3);
  alpha ~ normal(0, 3);
  n ~ multinomial(p);
}

generated quantities {
  vector[J] logit_psi;
  logit_psi = eta;
}

