data {
  int<lower=1> n;                 // PS sample size
  int<lower=1> n_np;              // NPS sample size
  int<lower=1> p;                 // number of coefficients

  matrix[n, p] X;                 // PS design matrix
  matrix[n_np, p] X_np;           // NPS design matrix

  vector<lower=0, upper=1>[n] y;        // PS outcomes
  vector<lower=0, upper=1>[n_np] y_np;   // NPS outcomes

  vector<lower=0>[n] w;           // PS weights (rescaled to sum to n)
  real<lower=0, upper=1> a;       // power prior exponent
}

parameters {
  vector[p] beta;
}

model {
  // Prior: t_3(0, 2.5)
  beta ~ student_t(3, 0, 2.5);

  // PS pseudolikelihood (vectorized)
  vector[n] eta = X * beta; // linear predictor
  target += dot_product(w, y .* eta - log1p_exp(eta));

  // NPS likelihood (power prior) (vectorized)
  vector[n_np] eta_np = X_np * beta;
  target += a * (y_np .* eta_np - log1p_exp(eta_np));
}
