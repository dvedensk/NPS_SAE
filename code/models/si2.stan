functions {
  // NOTE: sliced var (y_slice) comes FIRST, then start, end
  real bern_ll_slice(array[] int y_slice, int start, int end,
                     matrix X,
                     array[] int group,
                     vector alpha,
                     vector beta) {
    int M = end - start + 1;
    vector[M] eta;
    // matrix * vector → vector; add the random intercepts for this slice
    eta = X[start:end, ] * beta + alpha[group[start:end]];
    return bernoulli_logit_lpmf(y_slice | eta);
  }
}

// data {
//   int<lower=1> n;
//   int<lower=1> k;
//   matrix[n, k] X;
//   array[n] int<lower=0, upper=1> y;
//   int<lower=1> J;
//   array[n] int<lower=1, upper=J> group;
//   int<lower=1> grainsize;   // e.g., 250–1000
// }



// data {
//   // Training data
//   int<lower=1> n;
//   int<lower=1> k;
//   matrix[n, k] X;
//   array[n] int<lower=0, upper=1> y;
//   int<lower=1> J;
//   array[n] int<lower=1, upper=J> group;
//   int<lower=1> grainsize;
// 
//   // Ps population (probability-sample poststratification)
//   int<lower=0> Nps;
//   matrix[Nps, k] Xps;
//   array[Nps] int<lower=1, upper=J> group_ps;
// 
//   // P population (full population poststratification)
//   int<lower=0> Npop;
//   matrix[Npop, k] Xpop;
//   array[Npop] int<lower=1, upper=J> group_pop;
// }


data {
  // Training data
  int<lower=1> n;
  int<lower=1> k;
  matrix[n, k] X;
  array[n] int<lower=0, upper=1> y;
  int<lower=1> J;
  array[n] int<lower=1, upper=J> group;
  int<lower=1> grainsize;
  
  // ps
  int<lower=0> Nps;
  matrix[Nps, k] Xps;
  array[Nps] int<lower=1, upper=J> group_ps;
  vector<lower=0>[Nps] w_ps;


  // P
  int<lower=0> Npop;
  matrix[Npop, k] Xpop;
  array[Npop] int<lower=1, upper=J> group_pop;
  vector<lower=0>[Npop] w_pop;
}

parameters {
  vector[k] beta;
  real<lower=0> sigma_alpha;
  vector[J] alpha_raw;
}

transformed parameters {
  vector[J] alpha = sigma_alpha * alpha_raw;
}

model {
  // Priors
  beta ~ normal(0, 4);
  sigma_alpha ~ normal(0.2, 1);
  alpha_raw ~ normal(0, 1);

  // Parallelized likelihood
  target += reduce_sum(bern_ll_slice, y, grainsize, X, group, alpha, beta);
}


// using row
// generated quantities {
//   // PUMA-level weighted expectations
//   vector[J] mu_ps;   // weighted mean prob in Ps by PUMA
//   vector[J] mu_pop;  // weighted mean prob in P by PUMA
// 
//   // accumulators
//   vector[J] num_ps = rep_vector(0.0, J);
//   vector[J] den_ps = rep_vector(0.0, J);
//   vector[J] num_pop = rep_vector(0.0, J);
//   vector[J] den_pop = rep_vector(0.0, J);
// 
//   // Ps
//   for (i in 1:Nps) {
//     real eta = dot_product(row(Xps, i), beta) + alpha[group_ps[i]];
//     real p   = inv_logit(eta);
//     int j    = group_ps[i];
//     num_ps[j] += w_ps[i] * p;
//     den_ps[j] += w_ps[i];
//   }
//   for (j in 1:J)
//     mu_ps[j] = den_ps[j] > 0 ? num_ps[j] / den_ps[j] : negative_infinity();
// 
//   // P
//   for (i in 1:Npop) {
//     real eta = dot_product(row(Xpop, i), beta) + alpha[group_pop[i]];
//     real p   = inv_logit(eta);
//     int j    = group_pop[i];
//     num_pop[j] += w_pop[i] * p;
//     den_pop[j] += w_pop[i];
//   }
//   for (j in 1:J)
//     mu_pop[j] = den_pop[j] > 0 ? num_pop[j] / den_pop[j] : negative_infinity();
// }

// using cell
generated quantities {
  // ----- Ps (probability-sample poststratification; MRP-R) -----
  vector[J] mu_ps;                         // PUMA-level means (Ps-weighted)
  vector[J] num_ps = rep_vector(0.0, J);
  vector[J] den_ps = rep_vector(0.0, J);

  for (i in 1:Nps) {
    real eta = dot_product(row(Xps, i), beta) + alpha[group_ps[i]];
    real p   = inv_logit(eta);
    int j    = group_ps[i];
    num_ps[j] += w_ps[i] * p;              // weight by Ps weight (e.g., PWGTP)
    den_ps[j] += w_ps[i];
  }
  for (j in 1:J)
    mu_ps[j] = den_ps[j] > 0 ? num_ps[j] / den_ps[j] : not_a_number();

  // ----- Pop (full population poststratification; MRP-P) -----
  vector[J] mu_pop;                        // PUMA-level means (Pop-weighted)
  vector[J] num_pop = rep_vector(0.0, J);
  vector[J] den_pop = rep_vector(0.0, J);

  for (i in 1:Npop) {
    real eta = dot_product(row(Xpop, i), beta) + alpha[group_pop[i]];
    real p   = inv_logit(eta);
    int j    = group_pop[i];
    num_pop[j] += w_pop[i] * p;            // weight by cell size N_j (or 1s)
    den_pop[j] += w_pop[i];
  }
  for (j in 1:J)
    mu_pop[j] = den_pop[j] > 0 ? num_pop[j] / den_pop[j] : not_a_number();
}

