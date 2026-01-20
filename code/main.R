library(readr)
library(dplyr)
library(tidyverse)
library(sampling)
library(mvtnorm)
library(survey)
library(purrr)
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

nps_prior_mod <- cmdstan_model(
  file.path("code", "models", "nps_prior_logistic.stan"),
  cpp_options = list(stan_threads = TRUE)
)

Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())


set.seed(99)

# ============================================================
# CONFIGURATION: Response variable and sampling weights
# ============================================================
# Response variable to analyze (must exist in ACS_NPS_pop.csv)
response_var <- "PUBCOV" # Default: PUBCOV (binary). Can also use HICOV, WAGP, etc.
response_type <- "binary" # PUBCOV is 0/1 in processed population

# PS sampling weight configuration (for PUBCOV: use WAGP=0.05, PWGTP=-0.2)
PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)

# NPS sampling weight configuration (for PUBCOV: use PWGTP=0.3, AGEP=0.7)
NPS_weight_config <- list(PWGTP = 0.3, AGEP = 0.7)
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

alpha <- .05
mcmc_iter <- 1000
mcmc_burn <- 1000
n_chains <- 2
power_prior_a <- 0.5  # Power prior exponent for NPS Prior method (0-1)
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

# take samples
Nsim <- 1
prob_samples <- list()
nonprob_samples <- list()
results <- list()
summary_df_VSW <- list()
for (sim in 1:Nsim) {
  print(sim)

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
    sigma2_beta = 3
  )

  bulm_stan_out <- stan_bulm_mod$sample(
    data = bulm_stan_dat,
    chains = n_chains,
    parallel_chains = n_chains,
    iter_warmup = mcmc_burn,
    iter_sampling = mcmc_iter,
    threads_per_chain = 4
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
  weights_beta <- estimate_ipw(ps, nps, X_formula, "beta_reg")
  weights_uncond <- estimate_ipw(ps, nps, X_formula, "weighted")
  weights_beta <- c(weights_beta$ps_ipw, weights_beta$nps_ipw)
  weights_uncond <- c(weights_uncond$ps_ipw, weights_uncond$nps_ipw)
  scale_weights_beta <- weights_beta / sum(weights_beta) * length(y)
  scale_weights_uncond <- weights_uncond / sum(weights_uncond) * length(y)

  # 3.3 HT Direct Estimates (Horvitz-Thompson with IPW)
  beta_df <- data.frame(PUBCOV = y, PUMA = PUMAs, weights = weights_beta)
  uncond_df <- data.frame(PUBCOV = y, PUMA = PUMAs, weights = weights_uncond)
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
    X_formula = X_formula, Psi_formula = Psi_formula
  )
  bulm_beta_out$model <- "bulm_beta"

  bulm_uncond_out <- get_stan_summaries(
    y = y, X = X, Psi = Psi, weights = scale_weights_uncond,
    n_chains = n_chains, mcmc_burn = mcmc_burn,
    mcmc_iter = mcmc_iter, alpha = alpha,
    grouped_pop_df = acs_pop_grouped,
    X_formula = X_formula, Psi_formula = Psi_formula
  )
  bulm_uncond_out$model <- "bulm_uncond"

  # ============================================================
  # METHOD 4: MRP VARIANTS
  # ============================================================

  # 4.1 Basic MRP (MRP-R and MRP-P)
  mrp <- getMRP(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop
  )
  mrpr <- mrp$puma_summary_mrpr  # MRP-R (PS poststratification)
  mrpp <- mrp$puma_summary_mrpp  # MRP-P (Population poststratification)

  # 4.2 MRP with Integration (MRP-INT-R and MRP-INT-P)
  mrp1 <- getMRP_INT(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod = modINT,
    adjust = TRUE
  )
  mrpint_r <- mrp1$puma_summary_mrpr  # MRP-INT-R (PS poststratification)
  mrpint_p <- mrp1$puma_summary_mrpp  # MRP-INT-P (Population poststratification)

  # ============================================================
  # METHOD 5: VSW METHOD
  # ============================================================
  result_VSW <- vsw_out(ps[, !colnames(ps) %in% "weights"], nps, X_formula, response = response_var) 
  # This function returns a data frame with columns: "PUMA", "VSW_point_est", "ps_est", "nps_est", "pooled_results", "lower_CI", "upper_CI", "model"
  VSW_out <- result_VSW[,c("PUMA","VSW_point_est","lower_CI","upper_CI","model")]
  colnames(VSW_out) <- colnames(direst)

  # ============================================================
  # METHOD 6: NPS PRIOR METHOD (POWER PRIOR)
  # ============================================================
  # Build fixed-effects design matrices (X only, no spatial random effects)
  X_ps_fixed <- X_ps
  X_nps_fixed <- X_nps

  # Prepare Stan data
  nps_prior_data <- list(
    n = length(y_ps),
    n_np = length(y_nps),
    p = ncol(X_ps_fixed),
    X = X_ps_fixed,
    X_np = X_nps_fixed,
    y = y_ps,
    y_np = y_nps,
    w = ps_scale_weights,
    a = power_prior_a
  )

  # Sample from Stan model
  nps_prior_fit <- nps_prior_mod$sample(
    data = nps_prior_data,
    chains = n_chains,
    parallel_chains = n_chains,
    iter_warmup = mcmc_burn,
    iter_sampling = mcmc_iter,
    threads_per_chain = 4,
    refresh = 0  # Suppress progress messages
  )

  # Extract posterior draws of beta
  beta_draws <- posterior::as_draws_matrix(nps_prior_fit$draws("beta"))

  # Predict on population cells and aggregate by PUMA
  X_pop <- model.matrix(X_formula, data = acs_pop_grouped)
  n_cells <- nrow(X_pop)
  n_draws <- nrow(beta_draws)

  # Compute probabilities for each cell and draw
  logits <- X_pop %*% t(beta_draws)
  probs <- plogis(logits)

  # Generate predictions (binomial draws for each cell)
  preds <- matrix(
    rbinom(n_cells * n_draws, size = acs_pop_grouped$n, prob = c(probs)),
    nrow = n_cells, ncol = n_draws
  )

  # Aggregate by PUMA
  puma_preds <- acs_pop_grouped %>%
    select(PUMA, n) %>%
    mutate(pred_matrix = asplit(preds, 1)) %>%
    group_by(PUMA) %>%
    summarize(
      n_total = sum(n),
      pred_sums = list(colSums(do.call(rbind, pred_matrix))),
      .groups = "drop"
    ) %>%
    mutate(pred_probs = map(pred_sums, ~ .x / n_total)) %>%
    mutate(
      point_est = map_dbl(pred_probs, mean),
      lower_CI = map_dbl(pred_probs, ~ quantile(.x, alpha / 2)),
      upper_CI = map_dbl(pred_probs, ~ quantile(.x, 1 - alpha / 2))
    ) %>%
    select(PUMA, point_est, lower_CI, upper_CI) %>%
    mutate(model = "NPS Prior (Power)")

  nps_prior_out <- puma_preds

  # ============================================================
  # COMBINE ALL RESULTS
  # ============================================================
  results[[sim]] <- rbind(
    direst,              # METHOD 1
    bulm_ps_out,         # METHOD 2
    beta_HT,             # METHOD 3.3
    uncond_HT,           # METHOD 3.3
    bulm_beta_out,       # METHOD 3.4
    bulm_uncond_out,     # METHOD 3.4
    mrpr,                # METHOD 4.1
    mrpp,                # METHOD 4.1
    mrpint_r,            # METHOD 4.2
    mrpint_p,            # METHOD 4.2
    VSW_out,             # METHOD 5
    nps_prior_out        # METHOD 6
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

# VSW method doesn't have uncertainty quantification, so I return NA's for them. That'swhy I calculated it alone.

# Save results
save(
  list = c("prob_samples", "nonprob_samples", "results_df", "summary_df"),
  file = "data/ACS_NPS_simulation_results.RData"
)
