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
# load .stan if available (allows running non-Stan paths without failure)
stan_si2 <- file.path("code", "models", "si2.stan")
stan_mrp_int2 <- file.path("code", "models", "mrp_int2.stan")
stan_models_available <- file.exists(stan_si2) && file.exists(stan_mrp_int2)
if (!stan_models_available) {
  warning("Stan model files not found; skipping Stan compilation and Stan-based models.")
  mod <- modINT <- NULL
} else {
  mod <- cmdstan_model(
    stan_si2,
    cpp_options = list(stan_threads = TRUE)
  )

  modINT <- cmdstan_model(
    stan_mrp_int2,
    cpp_options = list(stan_threads = TRUE)
  )
  Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())
}


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

  # 5. Direct estimate on probability sample
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


  # 6. Fit unit-level model on probability sample

  bulm_out <- bulm_results(
    grouped_pop_df = acs_pop_grouped,
    alpha          = alpha,
    X              = X_ps,
    Psi            = Psi_ps,
    y              = y_ps,
    weights        = ps_scale_weights,
    sigma2_beta    = 1e4,
    iter           = 2000,
    burn           = 1000,
    X_formula      = X_formula,
    Psi_formula    = Psi_formula,
    summaries      = TRUE
  )
  bulm_out$model <- "bulm"

  # 7. Estimate IP weights for NPS and fit unit-level model
  #    on nonprobability sample with them (Tracy)
  ipw <- estimate_ipw(ps = ps, nps = nps, cov_formula = X_formula)


  bulm_ipw <- bulm_results(
    grouped_pop_df = acs_pop_grouped,
    alpha          = alpha,
    X              = X_nps,
    Psi            = Psi_nps,
    y              = y_nps,
    weights        = ipw,
    sigma2_beta    = 1e4,
    iter           = 2000,
    burn           = 1000,
    X_formula      = X_formula,
    Psi_formula    = Psi_formula,
    summaries      = TRUE
  )
  bulm_ipw$model <- "bulm_ipw"

  # 9. Combine results

  # 8. Fit NPS-informed prior model (Ethan)

  # 9. Fit MRP (Qianyu)



  mrp <- getMRP(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop
  )
  #  mrp_r
  mrpr <- mrp$puma_summary_mrpr

  #  mrp_p
  mrpp <- mrp$puma_summary_mrpp


  mrp1 <- getMRP_INT(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod = modINT
  )


  #  mrp_int
  mrpint <- mrp1$puma_summary_mrpp




  # 10. Fit VSW method (Qi)
  result_VSW <- vsw_out(ps[, !colnames(ps) %in% "weights"], nps, X_formula, response = response_var) 
  # This function returns a data frame with columns: "PUMA", "VSW_point_est", "ps_est", "nps_est", "pooled_results", "lower_CI", "upper_CI", "model"
  VSW_out <- result_VSW[,c("PUMA","VSW_point_est","lower_CI","upper_CI","model")]
  colnames(VSW_out) <- colnames(direst)
  # 11. Combine results

  results[[sim]] <- rbind(
    direst,
    bulm_out,
    bulm_ipw,
    VSW_out,
    mrpr,
    mrprp
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
