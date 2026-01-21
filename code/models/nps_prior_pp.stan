data {
  int<lower=1> n;                 // PS sample size
  int<lower=1> n_np;              // NPS sample size
  int<lower=1> p;                 // number of coefficients
  int<lower=1> J;                 // number of domains

  matrix[n, p] X;                 // PS design matrix
  matrix[n_np, p] X_np;           // NPS design matrix

  vector<lower=0, upper=1>[n] y;        // PS outcomes
  vector<lower=0, upper=1>[n_np] y_np;   // NPS outcomes

  vector<lower=0>[n] w;           // PS weights (rescaled to sum to n)
  real<lower=0, upper=1> a;       // power prior exponent

  array[n] int<lower=1, upper=J> domain;       // PS domain indices
  array[n_np] int<lower=1, upper=J> domain_np; // NPS domain indices
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
  // Prior: t_3(0, 2.5)
  beta ~ student_t(3, 0, 2.5);
  sigma_domain ~ normal(0, 1);
  z_domain ~ normal(0, 1);

  // PS pseudolikelihood (vectorized)
  vector[n] eta = X * beta + eta_domain[domain]; // linear predictor
  target += dot_product(w, y .* eta - log1p_exp(eta));

  // NPS likelihood (power prior) (vectorized)
  vector[n_np] eta_np = X_np * beta + eta_domain[domain_np];
  target += a * (y_np .* eta_np - log1p_exp(eta_np));
}
