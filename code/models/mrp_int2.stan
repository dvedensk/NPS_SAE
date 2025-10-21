functions {
  real bern_ll_slice(array[] int y_slice, int start, int end,
                     matrix X,
                     array[] int group,
                     vector alpha,
                     vector zeta,
                     array[] int psi_bin,
                     vector lp_psi,
                     vector beta,
                     real beta_psi) {
    int M = end - start + 1;
    vector[M] eta;
    // add group and psi-bin random intercepts + fixed logit(psi)
    eta =  X[start:end, ] * beta
         + beta_psi * lp_psi[start:end]
         + alpha[group[start:end]]
         + zeta[psi_bin[start:end]];
    return bernoulli_logit_lpmf(y_slice | eta);
  }
}

data {
  // --- Training data
  int<lower=1> n;
  int<lower=1> k;
  matrix[n, k] X;
  array[n] int<lower=0, upper=1> y;

  // group random intercept PUMA
  int<lower=1> J;
  array[n] int<lower=1, upper=J> group;

  // inclusion prob & bin for MRP-INT
  vector<lower=0, upper=1>[n] psi;          // estimated inclusion probabilities for training units
  int<lower=1> G;                           // number of psi bins
  array[n] int<lower=1, upper=G> psi_bin;   // bin index for each training unit

  // reduce_sum grainsize
  int<lower=1> grainsize;

  // --- Ps (probability-sample poststratification; MRP-R style)
  int<lower=0> Nps;
  matrix[Nps, k] Xps;
  array[Nps] int<lower=1, upper=J> group_ps;
  vector<lower=0, upper=1>[Nps] psi_ps;               // inclusion probs for Ps rows
  array[Nps] int<lower=1, upper=G> psi_bin_ps;        // corresponding psi-bin indices
  vector<lower=0>[Nps] w_ps;                          // e.g., survey weights

  // --- Pop (full population poststratification; MRP-P style)
  int<lower=0> Npop;
  matrix[Npop, k] Xpop;
  array[Npop] int<lower=1, upper=J> group_pop;
  vector<lower=0, upper=1>[Npop] psi_pop;             // inclusion probs for Pop rows
  array[Npop] int<lower=1, upper=G> psi_bin_pop;      // corresponding psi-bin indices
  vector<lower=0>[Npop] w_pop;                        // e.g., cell sizes N_j or 1s
}

transformed data {
  // precompute logit(psi) once (training only; Ps/Pop computed in GQ)
  vector[n] lp_psi = logit(psi);
}

parameters {
  vector[k] beta;
  real beta_psi;                     // fixed effect for logit(psi)

  real<lower=0> sigma_alpha;         // sd for group RE
  vector[J]     alpha_raw;

  real<lower=0> sigma_psi;           // sd for PUMA-bin RE
  vector[G]     zeta_raw;
}

transformed parameters {
  vector[J] alpha = sigma_alpha * alpha_raw;  // group RE
  vector[G] zeta  = sigma_psi  * zeta_raw;    //psi-bin RE
}

model {
  // Priors (weakly informative; tune as needed)
  beta       ~ normal(0, 4);
  beta_psi   ~ normal(0, 4);

  sigma_alpha ~ normal(0.2, 1);
  alpha_raw   ~ normal(0, 1);

  sigma_psi ~ normal(0.2, 1);
  zeta_raw  ~ normal(0, 1);

  // Parallelized likelihood (MRP-INT)
  target += reduce_sum(bern_ll_slice, y, grainsize,
                       X, group, alpha, zeta, psi_bin, lp_psi, beta, beta_psi);
}

generated quantities {
  // ----- Ps (probability-sample poststratification; MRP-R style summary) -----
  vector[J] mu_ps;                         // group-level weighted means
  vector[J] num_ps = rep_vector(0.0, J);
  vector[J] den_ps = rep_vector(0.0, J);

  // ----- Pop (full population poststratification; MRP-P style summary) -----
  vector[J] mu_pop;                        // group-level weighted means
  vector[J] num_pop = rep_vector(0.0, J);
  vector[J] den_pop = rep_vector(0.0, J);

  // Ps predictions (include psi fixed + psi-bin RE + group RE)
  for (i in 1:Nps) {
    real eta = dot_product(row(Xps, i), beta)
               + beta_psi * logit(psi_ps[i])
               + alpha[group_ps[i]]
               + zeta[psi_bin_ps[i]];
    real p   = inv_logit(eta);
    int j    = group_ps[i];
    num_ps[j] += w_ps[i] * p;
    den_ps[j] += w_ps[i];
  }
  for (j in 1:J)
    mu_ps[j] = den_ps[j] > 0 ? num_ps[j] / den_ps[j] : not_a_number();

  // Pop predictions (include psi fixed + psi-bin RE + group RE)
  for (i in 1:Npop) {
    real eta = dot_product(row(Xpop, i), beta)
               + beta_psi * logit(psi_pop[i])
               + alpha[group_pop[i]]
               + zeta[psi_bin_pop[i]];
    real p   = inv_logit(eta);
    int j    = group_pop[i];
    num_pop[j] += w_pop[i] * p;
    den_pop[j] += w_pop[i];
  }
  for (j in 1:J)
    mu_pop[j] = den_pop[j] > 0 ? num_pop[j] / den_pop[j] : not_a_number();
}
