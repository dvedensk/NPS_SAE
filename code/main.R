library(matrixStats)
library(tidyverse) # Loads readr, dplyr, purrr, ggplot2, etc.
library(sampling)
library(mvtnorm)
library(survey)
library(BayesLogit)
library(Matrix)
library(LaplacesDemon)
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
source(file.path("code", "utils.R")) # utils.R must define estimate_ipw()
source(file.path("code", "nps_prior.R")) # NPS prior helper functions
source(file.path("code", "models", "bulm.R"))
source(file.path("code", "models", "VSW.R"))

source(file.path("code", "models", "mrp_all.R"))

# Compile Stan models
mod <- cmdstan_model(
  file.path("code", "models", "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

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

Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())

# ============================================================
# CONFIGURATION: Set response & Data Defect Correlation
# ============================================================
# Response variable to analyze (must exist in ACS_NPS_pop.csv)
response_var <- "PUBCOV" # Default: PUBCOV (binary). Can also use HICOV, WAGP, etc.
response_type <- "binary" # PUBCOV is 0/1 in processed population

# PS sampling weight configuration (for PUBCOV: use WAGP=0.05, PWGTP=-0.2)
PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)

# NPS sampling weight configuration - POVPIP-based (avoids circularity with model covariates)
# Empirically grounded scenarios (Pew Research 2023 benchmarking study):
#   Extreme:   PWGTP=0.10, POVPIP=-1.52 → DDC=-0.089, ESS=6.6,   bias=-0.190 (worst opt-in error)
#   Typical:   PWGTP=0.10, POVPIP=-0.41 → DDC=-0.029, ESS=63.9,  bias=-0.061 (avg opt-in error)
#   Favorable: PWGTP=0.10, POVPIP=-0.22 → DDC=-0.014, ESS=254.6, bias=-0.030 (low-bias check)
NPS_weight_config <- list(PWGTP = 0.10, POVPIP = -1.52)  # Extreme (default)

# NPS prior configuration
# Choose from: "pp", "md", "mdl", "mdl10"
nps_prior_which <- c("pp")
scale_nps_prior_weight_covariate <- TRUE

# Other simulation parameters
Nsim <- 15
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

  # 1. Draw probability and nonprobability samples
  ps_sample <- get_strat_PS(pop_df = acs_pop, samp_frac = .005, weight_config = PS_weight_config)
  nps_sample <- get_NPS(pop_df = acs_pop, samp_frac = .05, weight_config = NPS_weight_config, internet_only = FALSE)

  # save the samples
  prob_samples[[sim]] <- ps_sample
  nonprob_samples[[sim]] <- nps_sample

  # check the PS weights
  cat("PS weight Check: sum(ps_weights) =", sum(ps_sample$weights), "vs N_pop =", N_pop, "\n")

  # Calculate and print sample diagnostics (DDC and ESS)
  diagnostics <- calculate_sample_diagnostics(
    pop_df = acs_pop,
    ps_sample = ps_sample,
    nps_sample = nps_sample,
    response_var = response_var,
    print_results = TRUE,
    sim_num = sim
  )

  # 2. Scale weights for pseudolikelihood models
  # Rescales weights to sum to sample size (for Bayesian pseudolikelihood)
  # Used by: BULM, NPS Prior (NOT used by design-based HT estimators)
  ps_scale_weights <- length(ps_sample$idx) * ps_sample$weights / sum(ps_sample$weights)
  nps_scale_weights <- length(nps_sample$idx) * nps_sample$weights / sum(nps_sample$weights)

  # 3. Extract PS and NPS data frames
  ps <- acs_pop[ps_sample$idx, ] %>%
    mutate(weights = ps_sample$weights)
  nps <- acs_pop[nps_sample$idx, ]

  # 4. Build design matrices for PS and NPS separately
  X_ps <- model.matrix(X_formula, data = ps)
  Psi_ps <- model.matrix(Psi_formula, data = ps)
  y_ps <- ps[[response_var]]

  X_nps <- model.matrix(X_formula, data = nps)
  Psi_nps <- model.matrix(Psi_formula, data = nps)
  y_nps <- nps[[response_var]]

  # ============================================================
  # METHOD 1: DIRECT ESTIMATES (PS ONLY)
  # ============================================================
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
      model = "direst"
    )

  # ============================================================
  # METHOD 2: STAN-BASED BULM (PS ONLY)
  # ============================================================
  bulm_stan_dat <- list(
    r = ncol(Psi_ps),
    nn = length(y_ps),
    p = ncol(X_ps),
    Y = y_ps,
    weights = ps_scale_weights,
    puma = apply(Psi_ps, 1, which.max),
    X = X_ps,
    sigma2_beta = 3 # Prior variance for fixed effects (weakly informative)
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
  bulm_ps_out$model <- "bulm_ps_only"

  # ============================================================
  # METHOD 3: IPW METHODS
  # ============================================================

  # 3.1 Data Preparation (Combine PS + NPS)
  X <- rbind(X_ps, X_nps)
  Psi <- rbind(Psi_ps, Psi_nps)
  y <- c(y_ps, y_nps)
  PUMAs <- c(ps$PUMA, nps$PUMA)

  # 3.2 IPW Weight Estimation (beta_reg & weighted methods)
  # TODO: Clarify difference between these methods with Dan (see OUTSTANDING_ISSUES.md)
  # - "beta_reg": Propensity score from beta regression?
  # - "weighted": Unconditional/marginal weighting?
  weights_beta <- estimate_ipw(ps, nps, X_formula, "beta_reg")
  weights_uncond <- estimate_ipw(ps, nps, X_formula, "weighted")
  weights_beta <- c(weights_beta$ps_ipw, weights_beta$nps_ipw)
  weights_uncond <- c(weights_uncond$ps_ipw, weights_uncond$nps_ipw)
  scale_weights_beta <- weights_beta / sum(weights_beta) * length(y)
  scale_weights_uncond <- weights_uncond / sum(weights_uncond) * length(y)

  # 3.3 HT Direct Estimates (Horvitz-Thompson with IPW)
  beta_df <- data.frame(response = y, PUMA = PUMAs, weights = weights_beta)
  names(beta_df)[1] <- response_var
  uncond_df <- data.frame(response = y, PUMA = PUMAs, weights = weights_uncond)
  names(uncond_df)[1] <- response_var
  beta_HT <- HT(beta_df, "beta_HT")
  uncond_HT <- HT(uncond_df, "uncond_HT")

  beta_HT$model <- "IPW HT (beta)"
  uncond_HT$model <- "IPW HT (uncond)"

  # 3.4 BULM with IPW Weights
  bulm_beta_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_beta,
    n_chains = n_chains, mcmc_burn = mcmc_burn,
    mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula,
    seed = curr_seed
  )
  bulm_beta_out$model <- "bulm_beta"

  bulm_uncond_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_uncond,
    n_chains = n_chains, mcmc_burn = mcmc_burn,
    mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula,
    seed = curr_seed
  )
  bulm_uncond_out$model <- "bulm_uncond"

  # ============================================================
  # METHOD 4: MRP VARIANTS
  # ============================================================
  # MRP Naming Convention:
  #   There is a MRP (plain) vs MRP-INT which integrates PS + NPS in joint model with selection correction
  #   Each of these MRP and MPR-INT has two different ways of poststratifying
  #     - R: Poststratifies on reweighted PS sample
  #     - P: Poststratifies on population frame (no bootstrap needed)
  # NOTE: adjust=TRUE --> PS weights scaled by (N_hat - n_np)/N_hat to account for NPS

  # 4.1 Basic MRP (MRP-R and MRP-P)
  mrp <- getMRP(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    bootstrap = TRUE, # Bootstrap needed for MRP-R uncertainty (not for MRP-P)
    L = 100,
    seed = curr_seed,
    n_chains = n_chains,
    stan_iter = mcmc_iter,
    stan_warmup = mcmc_burn,
    threads = mcmc_threads_per_chain
  )
  mrpr <- mrp$puma_summary_mrpr_bootstrap %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpp <- mrp$puma_summary_mrpp %>% select(PUMA, point_est, lower_CI, upper_CI, model)

  # 4.2 MRP with Integration (MRP-INT-R and MRP-INT-P)
  mrp1 <- getMRP_INT(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod = modINT,
    adjust = TRUE, # Adjust PS weights for population size (recommended)
    bootstrap = TRUE, # Bootstrap needed for MRP-INT-R uncertainty
    L = 100,
    seed = curr_seed,
    n_chains = n_chains,
    stan_iter = mcmc_iter,
    stan_warmup = mcmc_burn,
    threads = mcmc_threads_per_chain
  )
  mrpint_r <- mrp1$puma_summary_mrpr_bootstrap %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpint_p <- mrp1$puma_summary_mrpp %>% select(PUMA, point_est, lower_CI, upper_CI, model)

  
  # 4.3 MRP with Integration with PUBCOV (MRP-INT-R and MRP-INT-P)
  mrp_pubcov_out <- getMRP_new_PUBCOV(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod_int = mod_int,
    bootstrap = TRUE,
    L = 100,
    adjust = T,
    seed = curr_seed,
    n_chains = n_chains,
    stan_iter = mcmc_iter,
    stan_warmup = mcmc_burn,
    threads = mcmc_threads
  )
  mrpint_r_pubcov <- mrp_pubcov_out$puma_summary_mrpr_bootstrap %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  mrpint_p_pubcov <- mrp_pubcov_out$puma_summary_mrpp %>% select(PUMA, point_est, lower_CI, upper_CI, model)
  
  
  # ============================================================
  # METHOD 5: VSW METHOD
  # ============================================================
  # VSW: Variance-Scaled Weighting (see papers/aliste-et-al-2025-vsw-method.pdf)
  # NOTE: VSW does not produce interval estimates
  result_VSW <- vsw_out(ps[, !colnames(ps) %in% "weights"], nps, X_formula, response = response_var)
  VSW_out <- result_VSW[, c("PUMA", "VSW_point_est", "lower_CI", "upper_CI", "model")]
  colnames(VSW_out) <- colnames(direst)

  # ============================================================
  # METHOD 6: NPS PRIOR METHODS (MD/PP)
  # ============================================================
  domain_levels <- sort(unique(c(as.character(ps$PUMA), as.character(nps$PUMA))))
  domain_ps <- match(ps$PUMA, domain_levels)
  domain_nps <- match(nps$PUMA, domain_levels)
  PUMA_levels <- domain_levels

  power_prior_a <- NULL
  if ("pp" %in% nps_prior_which) {
    power_prior_a <- calculate_adaptive_power_prior_a(
      y_ps = y_ps,
      X_ps = X_ps,
      y_nps = y_nps,
      X_nps = X_nps,
      verbose = TRUE
    )$a
  }

  nps_prior_pl_res <- nps_prior_mcmc(
    md_mod = NULL,
    pp_mod = nps_prior_pp_mod,
    y = y_ps,
    X = X_ps,
    y_NP = y_nps,
    X_NP = X_nps,
    wts = ps_scale_weights,
    PUMA = factor(ps$PUMA, levels = domain_levels),
    PUMA_levels = PUMA_levels,
    raking_or_pl = "Pseudolikelihood",
    typeIerr = alpha,
    niter = mcmc_iter,
    warmup = mcmc_burn,
    chains = n_chains,
    seed = curr_seed,
    which_prior = nps_prior_which,
    a = power_prior_a,
    domain_ps = domain_ps,
    domain_nps = domain_nps,
    domain_levels = domain_levels,
    threads_per_chain = mcmc_threads_per_chain,
    parallel_chains = n_chains
  )

  # Raked weights for weight-as-covariate prior
  rake_vars <- c("AGEP_binned", "RAC1P", "SEX")
  pop_margins <- lapply(rake_vars, function(v) {
    as.data.frame(table(acs_pop[[v]])) %>%
      rename(!!v := Var1, Freq = Freq)
  })
  nps_design <- svydesign(ids = ~1, data = nps, weights = ~1)
  rake_formulas <- lapply(rake_vars, function(v) as.formula(paste0("~", v)))
  nps_raked_design <- rake(
    design = nps_design,
    sample.margins = rake_formulas,
    population.margins = pop_margins
  )
  nps_rake_weights <- weights(nps_raked_design)
  nps_rake_weights <- length(nps_rake_weights) * nps_rake_weights / sum(nps_rake_weights)

  nps_prior_rak_res <- nps_prior_mcmc(
    md_mod = NULL,
    pp_mod = nps_prior_pp_mod,
    y = y_ps,
    X = X_ps,
    y_NP = y_nps,
    X_NP = X_nps,
    wts = NULL,
    PUMA = factor(ps$PUMA, levels = domain_levels),
    PUMA_levels = PUMA_levels,
    raking_or_pl = "Raking",
    typeIerr = alpha,
    niter = mcmc_iter,
    warmup = mcmc_burn,
    chains = n_chains,
    seed = curr_seed,
    which_prior = nps_prior_which,
    a = power_prior_a,
    domain_ps = domain_ps,
    domain_nps = domain_nps,
    domain_levels = domain_levels,
    weight_covariate = list(ps = ps_scale_weights, nps = nps_rake_weights),
    scale_weight_covariate = scale_nps_prior_weight_covariate,
    threads_per_chain = mcmc_threads_per_chain,
    parallel_chains = n_chains
  )

  # ============================================================
  # COMBINE ALL RESULTS
  # ============================================================
  results[[sim]] <- rbind(
    direst, # METHOD 1
    bulm_ps_out, # METHOD 2
    beta_HT, # METHOD 3.3
    uncond_HT, # METHOD 3.3
    bulm_beta_out, # METHOD 3.4
    bulm_uncond_out, # METHOD 3.4
    mrpr, # METHOD 4.1
    mrpp, # METHOD 4.1
    mrpint_r, # METHOD 4.2
    mrpint_p, # METHOD 4.2
    mrpint_r_pubcov, # METHOD 4.3
    mrpint_p_pubcov, # METHOD 4.3
    VSW_out, # METHOD 5
    nps_prior_pl_res, # METHOD 6a
    nps_prior_rak_res # METHOD 6b
  )
}

# Aggregate across simulations
results_df <- results %>% list_rbind(names_to = "sim_num")

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

# VSW doesn't produce interval estimates, so Coverage and Interval Score will be NA

# Print results to console
cat(
  "\nSimulation complete. Summary by model (averaged over",
  Nsim, ifelse(Nsim == 1, "simulation", "simulations"), "):\n"
)
print(summary_df %>% arrange(MSE))

cat("\nExample PUMA-level estimates from the first simulation:\n")
print(
  results_df %>%
    filter(sim_num == min(as.numeric(sim_num))) %>%
    arrange(PUMA, model) %>%
    select(PUMA, model, point_est, lower_CI, upper_CI) %>%
    slice_head(n = 20)
)

# Save results
save(
  list = c("prob_samples", "nonprob_samples", "results_df", "summary_df"),
  file = "data/ACS_NPS_simulation_results.RData"
)
