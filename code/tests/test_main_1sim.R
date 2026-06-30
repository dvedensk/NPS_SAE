# Quick end-to-end test of main.R logic with Nsim=1
# Mirrors main.R exactly but overrides Nsim and uses fewer bootstrap samples.
# Purpose: confirm no errors in any code path before a full multi-sim run.

library(matrixStats)
library(tidyverse)
library(sampling)
library(mvtnorm)
library(survey)
library(BayesLogit)
library(Matrix)
library(LaplacesDemon)
library(cmdstanr)
tryCatch(cmdstanr::cmdstan_path(),
  error = function(e) stop("CmdStan not found.")
)

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R"))
source(file.path("code", "models", "nps_prior.R"))
source(file.path("code", "models", "bulm.R"))
source(file.path("code", "models", "VSW.R"))
source(file.path("code", "models", "ciginas.R"))
source(file.path("code", "models", "mrp_all.R"))

mod <- cmdstan_model(file.path("code", "models", "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)
modINT <- cmdstan_model(file.path("code", "models", "mrp_int2.stan"),
  cpp_options = list(stan_threads = TRUE)
)
stan_bulm_mod <- cmdstan_model(file.path("code", "models", "bulm.stan"),
  cpp_options = list(stan_threads = TRUE)
)
nps_prior_pp_mod <- cmdstan_model(file.path("code", "models", "nps_prior_pp.stan"),
  cpp_options = list(stan_threads = TRUE)
)

Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())

response_var <- "PUBCOV"
response_type <- "binary"
PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)
NPS_weight_config <- list(PWGTP = 0.10, POVPIP = -1.52)
nps_prior_which <- c("pp")
Nsim <- 1 # single replicate
alpha <- .05
mcmc_iter <- 1000
mcmc_burn <- 1000
n_chains <- 2
mcmc_threads_per_chain <- 4

acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv"), show_col_types = FALSE) %>%
  mutate(
    AGEP_binned = factor(AGEP_binned),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(response_true = mean(.data[[response_var]], na.rm = TRUE), .groups = "drop")

acs_pop_grouped <- acs_pop %>%
  group_by(PUMA, AGEP_binned, RAC1P, SEX) %>%
  tally()

X_formula <- as.formula("~ AGEP_binned + RAC1P + SEX")
Psi_formula <- as.formula("~ -1 + PUMA")
N_pop <- nrow(acs_pop)

prob_samples <- list()
nonprob_samples <- list()
results <- list()

for (sim in 1:Nsim) {
  curr_seed <- 99 + sim
  set.seed(curr_seed)
  cat("Starting simulation", sim, "with seed", curr_seed, "\n")

  ps_sample <- get_strat_PS(pop_df = acs_pop, samp_frac = .005, weight_config = PS_weight_config)
  nps_sample <- get_NPS(pop_df = acs_pop, samp_frac = .05, weight_config = NPS_weight_config, internet_only = FALSE)
  prob_samples[[sim]] <- ps_sample
  nonprob_samples[[sim]] <- nps_sample

  diagnostics <- calculate_sample_diagnostics(
    pop_df = acs_pop, ps_sample = ps_sample, nps_sample = nps_sample,
    response_var = response_var, print_results = TRUE, sim_num = sim
  )

  ps_scale_weights <- length(ps_sample$idx) * ps_sample$weights / sum(ps_sample$weights)
  nps_scale_weights <- length(nps_sample$idx) * nps_sample$weights / sum(nps_sample$weights)

  ps <- acs_pop[ps_sample$idx, ] %>% mutate(weights = ps_sample$weights)
  nps <- acs_pop[nps_sample$idx, ]

  X_ps <- model.matrix(X_formula, data = ps)
  Psi_ps <- model.matrix(Psi_formula, data = ps)
  y_ps <- ps[[response_var]]
  X_nps <- model.matrix(X_formula, data = nps)
  Psi_nps <- model.matrix(Psi_formula, data = nps)
  y_nps <- nps[[response_var]]

  # Method 1: Direct estimates
  samp.design <- svydesign(ids = ~1, weights = ~weights, data = ps)
  direst <- svyby(as.formula(paste0("~", response_var)), ~PUMA, samp.design, svymean,
    na.rm = TRUE, vartype = "se", keep.names = FALSE
  ) %>%
    arrange(PUMA) %>%
    transmute(PUMA,
      point_est = .data[[response_var]],
      lower_CI = point_est + qnorm(alpha / 2) * se,
      upper_CI = point_est + qnorm(1 - alpha / 2) * se,
      model = "direst"
    )

  # Method 2: BULM
  bulm_stan_dat <- list(
    r = ncol(Psi_ps), nn = length(y_ps), p = ncol(X_ps),
    Y = y_ps, weights = ps_scale_weights,
    puma = apply(Psi_ps, 1, which.max), X = X_ps, sigma2_beta = 3
  )
  bulm_stan_out <- stan_bulm_mod$sample(
    data = bulm_stan_dat, chains = n_chains,
    parallel_chains = n_chains, iter_warmup = mcmc_burn,
    iter_sampling = mcmc_iter,
    threads_per_chain = mcmc_threads_per_chain, seed = curr_seed
  )
  bulm_ps_out <- post_preds(
    grouped_pop_df = acs_pop_grouped,
    beta = bulm_stan_out$draws("beta"),
    eta = bulm_stan_out$draws("eta"),
    alpha = alpha, X_formula = X_formula,
    Psi_formula = Psi_formula, stan = TRUE
  )
  bulm_ps_out$model <- "bulm_ps_only"

  # Method 3: IPW
  X <- rbind(X_ps, X_nps)
  Psi <- rbind(Psi_ps, Psi_nps)
  y <- c(y_ps, y_nps)
  PUMAs <- c(ps$PUMA, nps$PUMA)
  weights_beta <- estimate_ipw(ps, nps, X_formula, "beta_reg")
  weights_uncond <- estimate_ipw(ps, nps, X_formula, "weighted")
  weights_beta <- c(weights_beta$ps_ipw, weights_beta$nps_ipw)
  weights_uncond <- c(weights_uncond$ps_ipw, weights_uncond$nps_ipw)
  scale_weights_beta <- weights_beta / sum(weights_beta) * length(y)
  scale_weights_uncond <- weights_uncond / sum(weights_uncond) * length(y)
  beta_df <- data.frame(response = y, PUMA = PUMAs, weights = weights_beta)
  names(beta_df)[1] <- response_var
  uncond_df <- data.frame(response = y, PUMA = PUMAs, weights = weights_uncond)
  names(uncond_df)[1] <- response_var
  beta_HT <- HT(beta_df, "beta_HT")
  beta_HT$model <- "IPW HT (beta)"
  uncond_HT <- HT(uncond_df, "uncond_HT")
  uncond_HT$model <- "IPW HT (uncond)"
  bulm_beta_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_beta,
    n_chains = n_chains, mcmc_burn = mcmc_burn, mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped, X_formula = X_formula, Psi_formula = Psi_formula, seed = curr_seed
  )
  bulm_beta_out$model <- "bulm_beta"
  bulm_uncond_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_uncond,
    n_chains = n_chains, mcmc_burn = mcmc_burn, mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped, X_formula = X_formula, Psi_formula = Psi_formula, seed = curr_seed
  )
  bulm_uncond_out$model <- "bulm_uncond"

  # Method 4: MRP
  mrp <- getMRP(
    MR = nps, ps = ps, acs_pop = acs_pop, mod = mod, bootstrap = TRUE, L = 100,
    seed = curr_seed, n_chains = n_chains, stan_iter = mcmc_iter, stan_warmup = mcmc_burn,
    threads = mcmc_threads_per_chain
  )
  mrpr <- mrp$puma_summary_mrpr_bootstrap %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpp <- mrp$puma_summary_mrpp %>% select(PUMA, point_est, lower_CI, upper_CI, model)

  mrp1 <- getMRP_INT(
    MR = nps, ps = ps, acs_pop = acs_pop, mod = modINT, adjust = TRUE,
    bootstrap = TRUE, L = 100, seed = curr_seed, n_chains = n_chains,
    stan_iter = mcmc_iter, stan_warmup = mcmc_burn, threads = mcmc_threads_per_chain
  )
  mrpint_r <- mrp1$puma_summary_mrpr_bootstrap %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpint_p <- mrp1$puma_summary_mrpp %>% select(PUMA, point_est, lower_CI, upper_CI, model)

  mrp_pubcov_out <- getMRP_INT(
    MR = nps, ps = ps, acs_pop = acs_pop, mod = modINT,
    include_response = TRUE, bootstrap = TRUE, L = 100, adjust = TRUE,
    seed = curr_seed, n_chains = n_chains, stan_iter = mcmc_iter,
    stan_warmup = mcmc_burn, threads = mcmc_threads_per_chain
  )
  mrpint_r_pubcov <- mrp_pubcov_out$puma_summary_mrpr_bootstrap %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpint_p_pubcov <- mrp_pubcov_out$puma_summary_mrpp %>% select(PUMA, point_est, lower_CI, upper_CI, model)

  # Method 5: VSW
  result_VSW <- vsw_out(ps[, !colnames(ps) %in% "weights"], nps, X_formula, response = response_var)
  VSW_out <- result_VSW[, c("PUMA", "VSW_point_est", "lower_CI", "upper_CI", "model")]
  colnames(VSW_out) <- colnames(direst)

  ciginas_res <- ciginas_out(
    ps = ps,
    nps = nps,
    X_formula = X_formula,
    response = response_var
  ) %>%
    select(PUMA, point_est, lower_CI, upper_CI, model)
  ciginas_res$model <- "Ciginas (simplified)"

  # Method 6: NPS Prior
  domain_levels <- sort(unique(c(as.character(ps$PUMA), as.character(nps$PUMA))))
  domain_ps <- match(ps$PUMA, domain_levels)
  domain_nps <- match(nps$PUMA, domain_levels)
  PUMA_levels <- domain_levels

  power_prior_a <- calculate_adaptive_power_prior_a(
    y_ps = y_ps, X_ps = X_ps,
    y_nps = y_nps, X_nps = X_nps, verbose = TRUE
  )$a

  nps_prior_pl_res <- nps_prior_mcmc(
    md_mod = NULL, pp_mod = nps_prior_pp_mod,
    y = y_ps, X = X_ps, y_NP = y_nps, X_NP = X_nps, wts = ps_scale_weights,
    PUMA = factor(ps$PUMA, levels = domain_levels), PUMA_levels = PUMA_levels,
    raking_or_pl = "Pseudolikelihood", typeIerr = alpha, niter = mcmc_iter, warmup = mcmc_burn,
    chains = n_chains, seed = curr_seed, which_prior = nps_prior_which, a = power_prior_a,
    domain_ps = domain_ps, domain_nps = domain_nps, domain_levels = domain_levels,
    threads_per_chain = mcmc_threads_per_chain, parallel_chains = n_chains
  )

  results[[sim]] <- rbind(
    direst, bulm_ps_out, beta_HT, uncond_HT, bulm_beta_out, bulm_uncond_out,
    mrpr, mrpp, mrpint_r, mrpint_p, mrpint_r_pubcov, mrpint_p_pubcov,
    VSW_out, ciginas_res, nps_prior_pl_res
  )
}

results_df <- results %>% list_rbind(names_to = "sim_num")

int_score <- function(alpha, y, lower, upper) {
  (upper - lower) +
    (2 / alpha) * (lower - y) * (y < lower) +
    (2 / alpha) * (y - upper) * (y > upper)
}

summary_df <- true_values %>%
  left_join(results_df, by = "PUMA") %>%
  group_by(model) %>%
  summarize(
    MSE = mean((response_true - point_est)^2),
    MAB = mean(abs(response_true - point_est)),
    Coverage = mean(between(response_true, lower_CI, upper_CI)),
    `Int. Score` = mean(int_score(alpha, response_true, lower_CI, upper_CI)),
    `Avg Width` = mean(upper_CI - lower_CI),
    .groups = "drop"
  )

cat("\n=== main.R 1-sim smoke test complete ===\n")
cat("Summary by model:\n")
print(summary_df %>% arrange(MSE), n = Inf)
