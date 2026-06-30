# Sensitivity analysis: sweep c in exp(-c * t²) for the exp-link power prior
# Run from project root: Rscript code/tests/sensitivity_c_sweep.R
#
# 5 sims × 3 c values = 15 MCMC runs, parallelized 3 at a time (12 cores)

library(matrixStats)
library(tidyverse)
library(sampling)
library(survey)
library(cmdstanr)

source("code/sampling_functions.R")
source("code/utils.R")
source("code/models/nps_prior.R")
source("code/models/bulm.R")

nps_prior_pp_mod <- cmdstan_model(
  "code/models/nps_prior_pp.stan",
  cpp_options = list(stan_threads = TRUE)
)
Sys.setenv(STAN_NUM_THREADS = 12)

acs_pop <- read_csv("data/ACS_NPS_pop.csv") %>%
  mutate(AGEP_binned = factor(AGEP_binned), RAC1P = factor(RAC1P),
         SEX = factor(SEX), PUMA = factor(PUMA))

true_values <- acs_pop %>% group_by(PUMA) %>%
  summarize(response_true = mean(PUBCOV), .groups = "drop")

acs_pop_grouped <- acs_pop %>% group_by(PUMA, AGEP_binned, RAC1P, SEX) %>% tally()
X_formula <- ~ AGEP_binned + RAC1P + SEX
Psi_formula <- ~ -1 + PUMA

Nsim <- 5
alpha <- 0.05
c_grid <- c(0.80, 0.90, 1.00)

# Pre-draw samples and compute t² for each sim (fast)
sim_data <- list()
for (sim in 1:Nsim) {
  curr_seed <- 99 + sim
  set.seed(curr_seed)
  cat("=== Drawing samples for sim", sim, "===\n")

  ps_sample <- get_strat_PS(pop_df = acs_pop, samp_frac = .005,
                             weight_config = list(WAGP = 0.05, PWGTP = -0.2))
  nps_sample <- get_NPS(pop_df = acs_pop, samp_frac = .05,
                         weight_config = list(PWGTP = 0.10, POVPIP = -1.52),
                         internet_only = FALSE)

  ps <- acs_pop[ps_sample$idx, ] %>% mutate(weights = ps_sample$weights)
  nps <- acs_pop[nps_sample$idx, ]
  ps_scale_weights <- length(ps_sample$idx) * ps_sample$weights / sum(ps_sample$weights)

  X_ps <- model.matrix(X_formula, data = ps)
  y_ps <- ps$PUBCOV
  X_nps <- model.matrix(X_formula, data = nps)
  y_nps <- nps$PUBCOV

  domain_levels <- sort(unique(c(as.character(ps$PUMA), as.character(nps$PUMA))))
  domain_ps <- match(ps$PUMA, domain_levels)
  domain_nps <- match(nps$PUMA, domain_levels)

  # Compute t² once (shared across c values)
  glm_ps  <- glm(y_ps  ~ X_ps  - 1, family = binomial())
  glm_nps <- glm(y_nps ~ X_nps - 1, family = binomial())
  diff_vec <- coef(glm_ps) - coef(glm_nps)
  V_sum <- vcov(glm_ps) + vcov(glm_nps)
  t2 <- as.numeric(t(diff_vec) %*% solve(V_sum) %*% diff_vec)

  sim_data[[sim]] <- list(
    y_ps = y_ps, X_ps = X_ps, y_nps = y_nps, X_nps = X_nps,
    ps_scale_weights = ps_scale_weights,
    domain_levels = domain_levels, domain_ps = domain_ps, domain_nps = domain_nps,
    ps = ps, t2 = t2, seed = curr_seed
  )
  cat("  t² =", round(t2, 4), "\n")
}

# Build all (sim, c) configs
configs <- expand.grid(sim = 1:Nsim, c_val = c_grid)
configs$a <- exp(-configs$c_val * sapply(configs$sim, function(s) sim_data[[s]]$t2))
cat("\n=== a values by (sim, c) ===\n")
print(configs %>% pivot_wider(names_from = c_val, values_from = a, names_prefix = "c="))

# Run MCMC in parallel (3 concurrent fits × 2 chains × 2 threads = 12 cores)
run_one <- function(i) {
  s <- configs$sim[i]
  c_val <- configs$c_val[i]
  d <- sim_data[[s]]
  a <- exp(-c_val * d$t2)

  res <- nps_prior_mcmc(
    d_mod = NULL, pp_mod = nps_prior_pp_mod,
    y = d$y_ps, X = d$X_ps, y_NP = d$y_nps, X_NP = d$X_nps,
    wts = d$ps_scale_weights,
    PUMA = factor(d$ps$PUMA, levels = d$domain_levels),
    PUMA_levels = d$domain_levels, raking_or_pl = "Pseudolikelihood",
    typeIerr = alpha, niter = 1000, warmup = 1000,
    chains = 2, seed = d$seed, which_prior = "pp", a = a,
    domain_ps = d$domain_ps, domain_nps = d$domain_nps,
    domain_levels = d$domain_levels,
    threads_per_chain = 2, parallel_chains = 2,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula
  )
  res$sim <- s
  res$c_val <- c_val
  res$a <- a
  res
}

cat("\n=== Running", nrow(configs), "MCMC fits (3 parallel) ===\n")
all_res <- parallel::mclapply(seq_len(nrow(configs)), run_one, mc.cores = 3)
results_df <- bind_rows(all_res)

# Compute performance metrics
summary_df <- true_values %>%
  left_join(results_df, by = "PUMA") %>%
  group_by(c_val, a) %>%
  summarize(
    MSE = mean((response_true - point_est)^2),
    MAB = mean(abs(response_true - point_est)),
    Coverage = mean(between(response_true, lower_CI, upper_CI)),
    .groups = "drop"
  )

# Aggregate: mean across sims for each c
c_summary <- summary_df %>%
  group_by(c_val) %>%
  summarize(
    mean_a = mean(a),
    MSE = mean(MSE),
    MAB = mean(MAB),
    Coverage = mean(Coverage),
    .groups = "drop"
  )

cat("\n\n========== C-VALUE SWEEP RESULTS (5 sims each) ==========\n")
print(c_summary)

saveRDS(list(results_df = results_df, summary_df = summary_df, c_summary = c_summary,
             configs = configs),
        "data/sensitivity_c_sweep.rds")
cat("\nResults saved to data/sensitivity_c_sweep.rds\n")
