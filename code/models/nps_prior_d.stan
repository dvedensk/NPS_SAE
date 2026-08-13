data {
  int<lower=1> n;                    // PS sample size
  int<lower=1> p;                    // number of predictors
  int<lower=1> J;                    // number of domains

  matrix[n, p] X;                    // PS design matrix

  array[n] int<lower=0, upper=1> y;  // PS outcomes

  vector<lower=0>[n] w;              // PS weights (rescaled to sum to n)

  vector[p] beta_prior_mean;         // beta hat NPS (will calculate in R)
  vector<lower=0>[p] beta_prior_sd;  // sigma_j

  array[n] int<lower=1, upper=J> domain; // PS domain indices
}

transformed data {
  vector[n] y_vec = to_vector(y);
}

parameters {
  vector[p] beta;
  vector[J] eta_domain;
  real<lower=0> sigma_domain;
}

model {
  // Prior
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  sigma_domain ~ cauchy(0, 5);
  eta_domain ~ normal(0, sigma_domain);

  // Pseudolikelihood (only PS) (vectorized)
  vector[n] lin_pred = X * beta + eta_domain[domain]; // linear predictor
  target += dot_product(w, y_vec .* lin_pred - log1p_exp(lin_pred));
}
