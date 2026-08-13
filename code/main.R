library(matrixStats)
library(tidyverse) # Loads readr, dplyr, purrr, ggplot2, etc.
library(sampling)
library(mvtnorm)
library(survey)
library(BayesLogit)
library(Matrix)
library(LaplacesDemon)
library(posterior)
library(rstanarm)

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr is required; install it via install.packages('cmdstanr').")
}
library(cmdstanr)
# Validate CmdStan is discoverable; cmdstanr 0.9.0 lacks cmdstan_available()
tryCatch(
  cmdstanr::cmdstan_path(),
  error = function(e) stop("CmdStan is not installed or toolchain is unavailable; run cmdstanr::install_cmdstan().")
)

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R")) 
source(file.path("code", "models", "nps_prior.R")) 
source(file.path("code", "models", "bulm.R"))
source(file.path("code", "models", "VSW.R"))
source(file.path("code", "models", "mrp_all.R"))
source(file.path("code", "mrp_helpers.R"))

# Compile Stan models
modINT <- cmdstan_model(
  file.path("code", "models", "mrp_int2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

stan_bulm_mod <- cmdstan_model(
  file.path("code", "models", "bulm.stan"),
  cpp_options = list(stan_threads = TRUE)
)

nps_prior_pp_mod <- cmdstan_model(
  file.path("code", "models", "nps_prior_pp.stan"),
  cpp_options = list(stan_threads = TRUE)
)

nps_prior_d_mod <- cmdstan_model(
  file.path("code", "models", "nps_prior_d.stan"),
  cpp_options = list(stan_threads = TRUE)
)

mrp_inclusion_mod <- cmdstan_model(
  file.path("code", "models", "mrp_int_propensity.stan"),
  cpp_options = list(stan_threads = TRUE)
)

mrp_inclusion_reff_mod <- cmdstan_model(
  file.path("code", "models", "mrp_int_propensity_random_effect.stan"),
  cpp_options = list(stan_threads = TRUE)
)

mrp_outcome_mod <- cmdstan_model(
  file.path("code", "models", "mrp_int_outcome.stan"),
  cpp_options = list(stan_threads = TRUE)
)

Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())

# ============================================================
# CONFIGURATION: Set response & Data Defect Correlation
# ============================================================
# Response variable to analyze (must exist in ACS_NPS_pop.csv)
response_var <- "PUBCOV" # Default: PUBCOV (binary). Can also use HICOV, WAGP, etc.
response_type <- "binary" # PUBCOV is 0/1 in processed population
setting <- "easy"

# PS sampling weight configuration (for PUBCOV: use WAGP=0.05, PWGTP=-0.2)
PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)

# NPS sampling weight configuration - POVPIP-based (avoids circularity with model covariates)
# Empirically grounded scenarios (Pew Research 2023 benchmarking study):
#   Extreme:   PWGTP=0.10, POVPIP=-1.52 → DDC=-0.089, ESS=6.6,   bias=-0.190 (worst opt-in error)
#   Typical:   PWGTP=0.10, POVPIP=-0.41 → DDC=-0.029, ESS=63.9,  bias=-0.061 (avg opt-in error)
#   Favorable: PWGTP=0.10, POVPIP=-0.22 → DDC=-0.014, ESS=254.6, bias=-0.030 (low-bias check)
if (setting == "easy") {
  NPS_weight_config <- list(PWGTP = 0.10, POVPIP = -0.22)
  save_dir <- file.path("results", "easy_setting_pubcov")
} else if (setting == "medium") {
  NPS_weight_config <- list(PWGTP = 0.10, POVPIP = -0.41)
  save_dir <- file.path("results", "medium_setting_pubcov")
} else if (setting == "hard") {
  NPS_weight_config <- list(PWGTP = 0.10, POVPIP = -1.52)
  save_dir <- file.path("results", "hard_setting_pubcov")
} else {
  stop("Improper setting")
}

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)


# Other simulation parameters
Nsim <- 100
alpha <- .05 # for interval estimates
mcmc_iter <- 1000
mcmc_burn <- 1000
n_chains <- 2
mcmc_threads_per_chain <- 4

# ============================================================
# LOAD POPULATION DATA
# ============================================================
# load population file
# Note: AGEP is continuous (for sampling weights), AGEP_binned is categorical (for modeling)
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv")) %>%
  mutate(
    AGEP_binned = factor(AGEP_binned),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# population values to compare against
true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    response_true = mean(.data[[response_var]], na.rm = TRUE),
    .groups = "drop"
  )

acs_pop_grouped <- acs_pop %>%
  group_by(PUMA, AGEP_binned, RAC1P, SEX) %>%
  tally()

X_formula <- as.formula("~ AGEP_binned + RAC1P + SEX")
X_formula_with_y <- update(X_formula, ~. +PUBCOV)
Psi_formula <- as.formula("~ -1 + PUMA")

N_pop <- nrow(acs_pop)
pop_mean <- mean(acs_pop[[response_var]], na.rm = TRUE)
cat(
  "\nPopulation summary:",
  "\n- N =", N_pop,
  "\n- Response mean (", response_var, ") =", round(pop_mean, 4),
  "\n- PS weights:", paste(names(PS_weight_config), PS_weight_config, sep = "=", collapse = ", "),
  "\n- NPS weights:", paste(names(NPS_weight_config), NPS_weight_config, sep = "=", collapse = ", "),
  "\n\n"
)

# ============================================================
# MAIN SIMULATION FOR-LOOP
# ============================================================
prob_samples <- list()
nonprob_samples <- list()
results <- list()
for (sim in 1:Nsim) {
  curr_seed <- 99 + sim
  set.seed(curr_seed) # Reproducible seed for each simulation
  cat("Starting simulation", sim, "with seed", curr_seed, "\n")

  # Draw probability and nonprobability samples
  ps_sample <- get_strat_PS(pop_df = acs_pop, samp_frac = .005,
                            weight_config = PS_weight_config)
  nps_sample <- get_NPS(pop_df = acs_pop, samp_frac = .05,
                        weight_config = NPS_weight_config, internet_only = FALSE)

  # save the samples
  prob_samples[[sim]] <- ps_sample
  nonprob_samples[[sim]] <- nps_sample

  # check the PS weights
  cat("PS weight Check: sum(ps_weights) =",
      sum(ps_sample$weights), "vs N_pop =", N_pop, "\n")

  # Calculate and print sample diagnostics (DDC and ESS)
  diagnostics <- calculate_sample_diagnostics(
    pop_df = acs_pop,
    ps_sample = ps_sample,
    nps_sample = nps_sample,
    response_var = response_var,
    print_results = TRUE,
    sim_num = sim
  )

  # Scale weights for pseudolikelihood models
  # Rescales weights to sum to sample size (for Bayesian pseudolikelihood)
  # Used by: BULM, NPS Prior (NOT used by design-based HT estimators)
  ps_scale_weights <- length(ps_sample$idx) * ps_sample$weights / sum(ps_sample$weights)
  nps_scale_weights <- length(nps_sample$idx) * nps_sample$weights / sum(nps_sample$weights)

  # Extract PS and NPS data frames
  ps <- acs_pop[ps_sample$idx, ] %>%
    mutate(weights = ps_sample$weights)
  nps <- acs_pop[nps_sample$idx, ]

  # Build design matrices for PS and NPS separately
  X_ps <- model.matrix(X_formula, data = ps)
  Psi_ps <- model.matrix(Psi_formula, data = ps)
  y_ps <- ps[[response_var]]

  X_nps <- model.matrix(X_formula, data = nps)
  Psi_nps <- model.matrix(Psi_formula, data = nps)
  y_nps <- nps[[response_var]]


  # MRP-INT-P with multinomial cell-level inclusion probabilities
  cell_vars <- c("PUMA", "AGEP_binned", "RAC1P", "SEX")
  
  mrp_cells <- make_mrp_inclusion_cells(
    pop_df = acs_pop,
    nps = nps,
    cell_vars = cell_vars
  )

  inc_fit <- fit_mrp_inclusion(
    cells = mrp_cells,
    X_formula = X_formula,
    puma_mode = "random",
    fixed_mod = mrp_inclusion_mod,
    random_mod = mrp_inclusion_reff_mod,
    chains = n_chains,
    nps=nps,
    cell_vars = cell_vars,
    response_var = response_var,
    model_name = "IPW-DE-CLIP",
    iter_warmup = mcmc_burn,
    iter_sampling = mcmc_iter,
    threads_per_chain = mcmc_threads_per_chain,
    seed = curr_seed
  )
  mrp_ipw_HT <- inc_fit$ipw_HT

  mrp_out_cells <- make_mrp_outcome_cells(
      cells = inc_fit$cells,
      nps = nps,
      response = response_var,
      cell_vars = cell_vars 
  )

  mrp_out_fit <- fit_mrp_outcome(
      cells = mrp_out_cells,
      X_formula = X_formula,
      mod = mrp_outcome_mod,
      chains = n_chains,
      iter_warmup = mcmc_burn,
      iter_sampling = mcmc_iter,
    threads_per_chain = mcmc_threads_per_chain,
      seed = curr_seed
  )

  mrp_est <- mrp_out_fit$summary_df
  mrp_est$model <- "MRP-INT-P"

  # Direct estimate on probability sample
  samp.design <- svydesign(ids = ~1, weights = ~weights, data = ps)

  # svyby extraction: use built-in se to avoid manual column juggling
  direst <- svyby(
    as.formula(paste0("~", response_var)),
    ~PUMA,
    samp.design,
    svymean,
    na.rm = TRUE,
    vartype = "se",
    keep.names = FALSE
  ) %>%
    arrange(PUMA) %>%
    transmute(
      PUMA,
      point_est = .data[[response_var]], # pulling the mean column from svyby
      lower_CI = point_est + qnorm(alpha / 2) * se,
      upper_CI = point_est + qnorm(1 - alpha / 2) * se,
      model = "Direct"
    )

  # Stan-based BULM (PS only)
  bulm_stan_dat <- list(
    r = ncol(Psi_ps),
    nn = length(y_ps),
    p = ncol(X_ps),
    Y = y_ps,
    weights = ps_scale_weights,
    puma = apply(Psi_ps, 1, which.max),
    X = X_ps,
    sigma2_beta = 9 
  )

  bulm_stan_out <- stan_bulm_mod$sample(
    data = bulm_stan_dat,
    chains = n_chains,
    parallel_chains = n_chains,
    iter_warmup = mcmc_burn,
    iter_sampling = mcmc_iter,
    threads_per_chain = mcmc_threads_per_chain,
    seed = curr_seed
  )

  bulm_ps_out <- post_preds(
    grouped_pop_df = acs_pop_grouped,
    beta = bulm_stan_out$draws("beta"),
    eta = bulm_stan_out$draws("eta"),
    alpha = alpha,
    X_formula = X_formula,
    Psi_formula = Psi_formula,
    stan = TRUE
  )
  bulm_ps_out$model <- "BULM-PS"

  # IPW Methods
  # Data Preparation (Combine PS + NPS)
  X <- rbind(X_ps, X_nps)
  Psi <- rbind(Psi_ps, Psi_nps)
  y <- c(y_ps, y_nps)
  PUMAs <- c(ps$PUMA, nps$PUMA)

  # IPW Weight Estimation 
  weights_uncond <- estimate_ipw(ps, nps, X_formula, "weighted")
  weights_uncond <- c(weights_uncond$ps_ipw, weights_uncond$nps_ipw)
  scale_weights_uncond <- weights_uncond / sum(weights_uncond) * length(y)

  weights_uncond_y <- estimate_ipw(ps, nps, X_formula_with_y, "weighted")
  weights_uncond_y <- c(weights_uncond_y$ps_ipw, weights_uncond_y$nps_ipw)
  scale_weights_uncond_y <- weights_uncond_y / sum(weights_uncond_y) * length(y)

  # IPW Direct Estimates
  uncond_df <- data.frame(response = y, PUMA = PUMAs, weights = weights_uncond)
  names(uncond_df)[1] <- response_var
  uncond_HT <- HT(uncond_df, "uncond_HT")

  uncond_y_df <- data.frame(response = y, PUMA = PUMAs, weights = weights_uncond_y)
  names(uncond_y_df)[1] <- response_var
  uncond_y_HT <- HT(uncond_y_df, "uncond_y_HT")

  uncond_HT$model <- "IPW-DE"
  uncond_y_HT$model <- "IPW-DE+Y"

  # BULM with IPW Weights
  bulm_uncond_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_uncond,
    n_chains = n_chains, mcmc_burn = mcmc_burn,
    mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula,
    seed = curr_seed
  )
  bulm_uncond_out$model <- "IPW-BULM"

  bulm_uncond_y_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_uncond_y,
    n_chains = n_chains, mcmc_burn = mcmc_burn,
    mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula,
    seed = curr_seed
  )
  bulm_uncond_y_out$model <- "IPW-BULM+Y"

  # MRP with Integration and combined-sample weights
  mrp1 <- getMRP_INT(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod = modINT,
    adjust = TRUE, # Adjust PS weights for population size (recommended)
    bootstrap = FALSE, # Bootstrap needed for MRP-INT-R uncertainty
    L = 100,
    seed = curr_seed,
    n_chains = n_chains,
    stan_iter = mcmc_iter,
    stan_warmup = mcmc_burn,
    threads = mcmc_threads_per_chain
  )
  mrpint_p <- mrp1$puma_summary_mrpp %>%
      select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpint_p$model <- "MRP-INT-P-IPW"

  # VSW method (does not produce interval estimates)
  result_VSW <- vsw_out(ps[, !colnames(ps) %in% "weights"], nps, X_formula, response = response_var)
  VSW_out <- result_VSW[, c("PUMA", "VSW_point_est", "lower_CI", "upper_CI", "model")]
  colnames(VSW_out) <- colnames(direst)

  # NPS prior methods (md/pp)
  domain_levels <- sort(unique(c(as.character(ps$PUMA), as.character(nps$PUMA))))
  domain_ps <- match(ps$PUMA, domain_levels)
  domain_nps <- match(nps$PUMA, domain_levels)
  PUMA_levels <- domain_levels

  # Compute both a values: p-value link and exp-link
  pp_pval <- calculate_adaptive_power_prior_a(y_ps, X_ps, y_nps, X_nps, verbose = TRUE)
  pp_exp  <- calculate_exp_link_a(y_ps, X_ps, y_nps, X_nps, null_target = 0.9, verbose = TRUE)

  # Shared args for nps_prior_mcmc
  nps_prior_shared <- list(
    d_mod = nps_prior_d_mod, pp_mod = nps_prior_pp_mod,
    y = y_ps, X = X_ps, y_NP = y_nps, X_NP = X_nps,
    wts = ps_scale_weights,
    PUMA = factor(ps$PUMA, levels = domain_levels),
    PUMA_levels = PUMA_levels, raking_or_pl = "Pseudolikelihood",
    typeIerr = alpha, niter = mcmc_iter, warmup = mcmc_burn,
    chains = n_chains, seed = curr_seed, which_prior = "pp",
    domain_ps = domain_ps, domain_nps = domain_nps,
    domain_levels = domain_levels,
    threads_per_chain = mcmc_threads_per_chain, parallel_chains = n_chains,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula
  )

  # Power prior with p-value link
  nps_prior_pval_res <- do.call(nps_prior_mcmc, c(nps_prior_shared, list(a = pp_pval$a)))
  nps_prior_pval_res$model <- "NIP (p-value)"

  # Power prior with exp-link
  nps_prior_exp_res <- do.call(nps_prior_mcmc, c(nps_prior_shared, list(a = pp_exp$a)))
  nps_prior_exp_res$model <- "NIP (exp link)"

  # Combine all results
  results[[sim]] <- rbind(
    direst, 
    bulm_ps_out, 
    uncond_HT, 
    bulm_uncond_out, 
    uncond_y_HT, 
    bulm_uncond_y_out, 
    mrp_ipw_HT,
    mrp_est,
    mrpint_p,
    VSW_out, 
    nps_prior_pval_res, 
    nps_prior_exp_res 
  )

  results[[sim]]$sim_num <- sim
  sim_results <- results[[sim]]

  save(
    sim_results,
    file = file.path(save_dir, paste0("intermediate_result_", sim, ".RData"))
  )
}

# Aggregate across simulations
results_df <- results %>% list_rbind()

# Summaries for each method
summary_df <- true_values %>%
  left_join(results_df, by = "PUMA") %>%
  group_by(model) %>%
  summarize(
    MSE = mean((response_true - point_est)^2),
    MAB = mean(abs(response_true - point_est)),
    Coverage = mean(between(response_true, lower_CI, upper_CI)),
    `Int. Score` = mean(int_score(alpha, response_true, lower_CI, upper_CI))
  )

save(
  list = c("prob_samples", "nonprob_samples", "results_df", "summary_df"),
  file = file.path(save_dir, "ACS_NPS_simulation_results.RData")
)
