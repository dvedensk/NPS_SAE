// LEGACY: older NPS prior model (logistic, no domain effects).
// Current implementation uses nps_prior_pp.stan; kept for reference only.

data {
  int<lower=1> n;                 // PS sample size
  int<lower=1> n_np;              // NPS sample size
  int<lower=1> p;                 // number of coefficients

  matrix[n, p] X;                 // PS design matrix
  matrix[n_np, p] X_np;           // NPS design matrix

  array[n] int<lower=0, upper=1> y;        // PS outcomes
  array[n_np] int<lower=0, upper=1> y_np;  // NPS outcomes

  vector<lower=0>[n] w;           // PS weights (rescaled to sum to n)
  real<lower=0, upper=1> a;       // power prior exponent
}

parameters {
  vector[p] beta;
}

model {
  // Prior: t_3(0, 2.5)
  beta ~ student_t(3, 0, 2.5);

  // PS pseudolikelihood
  for (i in 1:n) {
    target += w[i] * bernoulli_logit_lpmf(y[i] | X[i] * beta);
  }

  // NPS likelihood (power prior)
  for (l in 1:n_np) {
    target += a * bernoulli_logit_lpmf(y_np[l] | X_np[l] * beta);
  }
}
