data {
  int<lower=1> n;                    // PS sample size
  int<lower=1> p;                    // number of predictors

  matrix[n, p] X;                    // PS design matrix

  vector<lower=0, upper=1>[n] y;  // PS outcomes

  vector<lower=0>[n] w;              // PS weights (rescaled to sum to n)

  vector[p] beta_prior_mean;         // beta hat NPS (will calculate in R)
  vector<lower=0>[p] beta_prior_sd;  // sigma_j
}

parameters {
  vector[p] beta;
}

model {
  // Prior
  beta ~ normal(beta_prior_mean, beta_prior_sd);

  // Pseudolikelihood (only PS) (vectorized)
  vector[n] eta = X * beta; // linear predictor
  target += dot_product(w, y .* eta - log1p_exp(eta));
}
