data {
  int<lower=1> n;                 // PS sample size
  int<lower=1> n_np;              // NPS sample size
  int<lower=1> p;                 // number of coefficients
  int<lower=1> J;                 // number of domains

  matrix[n, p] X;                 // PS design matrix
  matrix[n_np, p] X_np;           // NPS design matrix

  array[n] int<lower=0, upper=1> y;        // PS outcomes
  array[n_np] int<lower=0, upper=1> y_np;   // NPS outcomes

  vector<lower=0>[n] w;           // PS weights (rescaled to sum to n)
  real<lower=0, upper=1> a;       // power prior exponent

  array[n] int<lower=1, upper=J> domain;       // PS domain indices
  array[n_np] int<lower=1, upper=J> domain_np; // NPS domain indices
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
  beta ~ normal(0, 3);
  sigma_domain ~ cauchy(0, 5);
  eta_domain ~ normal(0, sigma_domain);

  // PS pseudolikelihood 
  vector[n] eta = X * beta + eta_domain[domain];
  target += dot_product(w, y_vec .* eta - log1p_exp(eta));

  // NPS likelihood (power prior)
  if(a > 0) {
    target += a * bernoulli_logit_glm_lpmf(y_np | X_np, eta_domain[domain_np], beta);
  }
}
