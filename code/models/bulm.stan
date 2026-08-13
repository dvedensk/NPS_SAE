data {
  int<lower=1> r; // Number of counties
  int<lower=0> nn; // Number of observations
  int<lower=1> p; // Dimension of fixed effect coefficients
  array[nn] int<lower=0, upper=1> Y; // Responses
  vector[nn] weights;  // weights
  array[nn] int<lower=1, upper=r> puma; // County number of each observation
  matrix[nn, p] X; // Fixed effect model matrix
  real<lower=0> sigma2_beta; // fixed effect variance
}
parameters {
  vector[r] eta; // random effects
  real<lower=0> sigma_eta; // standard deviation of random effects
  vector[p] beta; // fixed effects
}
model {
  vector[nn] linPred;
  eta ~ normal(0, sigma_eta);
  sigma_eta ~ cauchy(0, 5);
  linPred = X*beta + eta[puma];
  for(i in 1:nn){
    target += weights[i] * bernoulli_logit_lpmf(Y[i] | linPred[i]); 
  }
  for(b in 1:p){
    beta[b] ~ normal(0,sqrt(sigma2_beta));
  }
}
