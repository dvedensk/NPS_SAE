data {
  int<lower=1> n;                    // PS sample size
  int<lower=1> p;                    // number of predictors
  int<lower=1> J;                    // number of domains

  matrix[n, p] X;                    // PS design matrix

  vector<lower=0, upper=1>[n] y;  // PS outcomes

  vector<lower=0>[n] w;              // PS weights (rescaled to sum to n)

  vector[p] beta_prior_mean;         // beta hat NPS (will calculate in R)
  vector<lower=0>[p] beta_prior_sd;  // sigma_j

  array[n] int<lower=1, upper=J> domain; // PS domain indices
}

parameters {
  vector[p] beta;
  vector[J] z_domain;
  real<lower=0> sigma_domain;
}

transformed parameters {
  vector[J] eta_domain;
  eta_domain = sigma_domain * z_domain; // non-centered for better efficiency
}

model {
  // Prior
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  sigma_domain ~ cauchy(0, 5);
  z_domain ~ normal(0, 1);

  // Pseudolikelihood (only PS) (vectorized)
  vector[n] eta = X * beta + eta_domain[domain]; // linear predictor
  target += dot_product(w, y .* eta - log1p_exp(eta));
}
