library(dplyr)
library(readr)
library(sampling)
library(survey)

source(file.path("code", "sampling_functions.R"))

# ============================================================
# CONFIGURATION: Response variable and sampling weights
# ============================================================
# Response variable to analyze (must exist in ACS_NPS_pop.csv)
RESPONSE_VAR <- "PUBCOV" # Default: PUBCOV (binary). Can also use HICOV, WAGP, etc.

# PS sampling weight configuration (for PUBCOV: use WAGP=0.05, PWGTP=-0.2)
PS_WEIGHT_CONFIG <- list(WAGP = 0.05, PWGTP = -0.2)

# NPS sampling weight configuration (for PUBCOV: use PWGTP=0.3, AGEP=0.7)
NPS_WEIGHT_CONFIG <- list(PWGTP = 0.3, AGEP = 0.7)
# ============================================================

# Load population data
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv"))

# Population info
total_pumas <- length(unique(acs_pop$PUMA))

# Sample fractions (NPS is 10x larger than PS)
PS_SAMP_FRAC <- 0.005 # 0.5% → ~7,378 individuals
NPS_SAMP_FRAC <- 0.05 # 5% → ~73,781 individuals (10x PS)

# Number of replications
N_REPS <- 20

# Optional: Check for zero-variance PUMAs (adds runtime)
CHECK_ZERO_VAR <- TRUE

cat("Testing correlations with", RESPONSE_VAR, ":\n")
cat("PS weight_config:", paste(names(PS_WEIGHT_CONFIG), "=", PS_WEIGHT_CONFIG, collapse = ", "), "\n")
cat("NPS weight_config:", paste(names(NPS_WEIGHT_CONFIG), "=", NPS_WEIGHT_CONFIG, collapse = ", "), "\n")
cat("PS sample fraction:", PS_SAMP_FRAC, "\n")
cat("NPS sample fraction:", NPS_SAMP_FRAC, "\n")
cat("Check zero-variance PUMAs:", CHECK_ZERO_VAR, "\n")
cat("Running", N_REPS, "replications...\n\n")

# Storage for results
ps_correlations <- numeric(N_REPS)
nps_correlations <- numeric(N_REPS)
ps_ess <- numeric(N_REPS) # Effective sample size
nps_ess <- numeric(N_REPS)
ps_puma_coverage <- numeric(N_REPS) # Number of PUMAs covered
nps_puma_coverage <- numeric(N_REPS)
ps_response_mean <- numeric(N_REPS) # Mean response variable in sample
nps_response_mean <- numeric(N_REPS)
if (CHECK_ZERO_VAR) {
  ps_zero_var_counts <- numeric(N_REPS)
}

for (i in 1:N_REPS) {
  set.seed(i)

  # Draw PS sample
  ps <- get_strat_PS(
    pop_df = acs_pop,
    samp_frac = PS_SAMP_FRAC,
    weight_config = PS_WEIGHT_CONFIG
  )

  # Draw NPS sample (using defaults: internet_only = FALSE for maximum bias)
  nps <- get_NPS(
    pop_df = acs_pop,
    samp_frac = NPS_SAMP_FRAC,
    weight_config = NPS_WEIGHT_CONFIG
    # internet_only defaults to FALSE
  )

  # Calculate DDC and ESS using helper function
  diagnostics <- calculate_sample_diagnostics(
    pop_df = acs_pop,
    ps_sample = ps,
    nps_sample = nps,
    response_var = RESPONSE_VAR,
    print_results = FALSE
  )

  # Store results
  ps_correlations[i] <- diagnostics$ps_ddc
  nps_correlations[i] <- diagnostics$nps_ddc
  ps_ess[i] <- diagnostics$ps_ess
  nps_ess[i] <- diagnostics$nps_ess

  # Calculate PUMA coverage (number of unique PUMAs in each sample)
  ps_puma_coverage[i] <- length(unique(acs_pop$PUMA[ps$idx]))
  nps_puma_coverage[i] <- length(unique(acs_pop$PUMA[nps$idx]))

  # Calculate mean response variable in each sample
  ps_response_mean[i] <- mean(acs_pop[[RESPONSE_VAR]][ps$idx])
  nps_response_mean[i] <- mean(acs_pop[[RESPONSE_VAR]][nps$idx])

  # Check for zero-variance direct estimates by PUMA for the PS only (optional)
  if (CHECK_ZERO_VAR) {
    # PS sample
    ps_sample <- acs_pop[ps$idx, ]
    ps_sample$design_weight <- ps$weights
    ps_design <- svydesign(ids = ~1, weights = ~design_weight, data = ps_sample)
    ps_formula <- as.formula(paste0("~", RESPONSE_VAR))
    ps_direst <- svyby(ps_formula, ~PUMA, ps_design, svymean, vartype = "var", na.rm = TRUE)
    # Robustly extract variance column (handles both "var" and "var.RESPONSE_VAR" naming)
    var_cols <- names(ps_direst)[grepl("^var", names(ps_direst))]
    ps_zero_var_counts[i] <- sum(ps_direst[[var_cols[1]]] == 0, na.rm = TRUE)
  }

  if (i %% 10 == 0) {
    cat("Completed", i, "replications...\n")
  }
}

# Calculate population mean response variable
pop_response_mean <- mean(acs_pop[[RESPONSE_VAR]])

cat("\n========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n")

cat("\n--- Population ---\n")
cat(sprintf("%s: %.1f%% | PUMAs: %d\n", RESPONSE_VAR, 100 * pop_response_mean, total_pumas))

cat("\n--- PS Sample ---\n")
cat(sprintf("DDC: %.4f (SD: %.4f)\n", mean(ps_correlations), sd(ps_correlations)))
cat(sprintf(
  "ESS: %d (median: %d, range: %d-%d)\n",
  round(mean(ps_ess)), round(median(ps_ess)), round(min(ps_ess)), round(max(ps_ess))
))
cat(sprintf(
  "PUMA Coverage: %.1f/%.0f (%.1f%%)\n",
  mean(ps_puma_coverage), total_pumas, 100 * mean(ps_puma_coverage) / total_pumas
))
cat(sprintf(
  "Sample %s: %.1f%% (bias: %+.2f pp)\n",
  RESPONSE_VAR, 100 * mean(ps_response_mean), 100 * (mean(ps_response_mean) - pop_response_mean)
))
if (CHECK_ZERO_VAR) {
  cat(sprintf(
    "Zero-var PUMAs: %.1f (%.1f%%)\n",
    mean(ps_zero_var_counts), 100 * mean(ps_zero_var_counts) / total_pumas
  ))
}

cat("\n--- NPS Sample ---\n")
cat(sprintf(
  "DDC: %.3f (SD: %.4f, range: %.3f to %.3f)\n",
  mean(nps_correlations), sd(nps_correlations), min(nps_correlations), max(nps_correlations)
))
cat(sprintf(
  "ESS: %.1f (median: %.1f, range: %.1f-%.1f)\n",
  mean(nps_ess), median(nps_ess), min(nps_ess), max(nps_ess)
))
cat(sprintf(
  "PUMA Coverage: %.1f/%.0f (%.1f%%)\n",
  mean(nps_puma_coverage), total_pumas, 100 * mean(nps_puma_coverage) / total_pumas
))
cat(sprintf(
  "Sample %s: %.1f%% (bias: %+.2f pp)\n",
  RESPONSE_VAR, 100 * mean(nps_response_mean), 100 * (mean(nps_response_mean) - pop_response_mean)
))

cat("\n========================================\n")
