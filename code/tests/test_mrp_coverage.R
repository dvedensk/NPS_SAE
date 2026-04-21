#!/usr/bin/env Rscript
# Comprehensive MRP coverage test: Compare baseline vs bootstrap for MRP and MRP-INT

library(readr)
library(dplyr)
library(cmdstanr)
library(here)
library(survey)

setwd(here::here())
source("code/tests/create_test_sample.R")
source("code/models/mrp_all.R")

# Load population
acs_pop <- read_csv("data/ACS_NPS_pop.csv") %>%
  mutate(
    AGEP_binned = factor(AGEP_binned),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

cat("\n=== MRP Coverage Comparison: Baseline vs Bootstrap ===\n")
cat("Testing MRP and MRP-INT (with and without include_response)\n\n")

# Compile Stan models
cat("Compiling Stan models...\n")
mod <- cmdstan_model(
  here("code/models", "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)
cat("✓ MRP model compiled\n")

mod_int <- cmdstan_model(
  here("code/models", "mrp_int2.stan"),
  cpp_options = list(stan_threads = TRUE)
)
cat("✓ MRP-INT model compiled\n\n")

# Calculate true values
response_var <- "PUBCOV"
alpha <- 0.05

true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    response_true = mean(.data[[response_var]], na.rm = TRUE),
    .groups = "drop"
  )

# Direct estimates (for comparison)
cat("Computing direct estimates...\n")
samp.design <- svydesign(ids = ~1, weights = ~weights, data = ps_df)
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
    point_est = .data[[response_var]],
    lower_CI = point_est + qnorm(alpha / 2) * se,
    upper_CI = point_est + qnorm(1 - alpha / 2) * se,
    model = "direst"
  )

# ============================================================
# TEST 1: MRP (baseline, no bootstrap)
# ============================================================
cat("\n[1/4] Running MRP (baseline, no bootstrap)...\n")
start_time <- Sys.time()
mrp_baseline <- getMRP(
  MR = nps_df,
  ps = ps_df,
  acs_pop = acs_pop,
  mod = mod,
  bootstrap = FALSE
)
elapsed <- Sys.time() - start_time
cat("✓ Completed in", round(elapsed, 2), attr(elapsed, "units"), "\n")

mrp_r_baseline <- mrp_baseline$puma_summary_mrpr
mrp_p_baseline <- mrp_baseline$puma_summary_mrpp

# ============================================================
# TEST 2: MRP (with bootstrap, L=100)
# ============================================================
cat("\n[2/4] Running MRP (with bootstrap, L=100)...\n")
start_time <- Sys.time()
mrp_boot <- getMRP(
  MR = nps_df,
  ps = ps_df,
  acs_pop = acs_pop,
  mod = mod,
  bootstrap = TRUE,
  L = 100
)
elapsed <- Sys.time() - start_time
cat("✓ Completed in", round(elapsed, 2), attr(elapsed, "units"), "\n")

mrp_r_boot <- mrp_boot$puma_summary_mrpr_bootstrap

# ============================================================
# TEST 3: MRP-INT (baseline, no bootstrap)
# ============================================================
cat("\n[3/4] Running MRP-INT (baseline, no bootstrap)...\n")
start_time <- Sys.time()
mrp_int_baseline <- getMRP_INT(
  MR = nps_df,
  ps = ps_df,
  acs_pop = acs_pop,
  mod = mod_int,
  adjust = TRUE,
  bootstrap = FALSE
)
elapsed <- Sys.time() - start_time
cat("✓ Completed in", round(elapsed, 2), attr(elapsed, "units"), "\n")

mrp_int_r_baseline <- mrp_int_baseline$puma_summary_mrpr
mrp_int_p_baseline <- mrp_int_baseline$puma_summary_mrpp

# ============================================================
# TEST 4: MRP-INT (with bootstrap, L=100)
# ============================================================
cat("\n[4/4] Running MRP-INT (with bootstrap, L=100)...\n")
start_time <- Sys.time()
mrp_int_boot <- getMRP_INT(
  MR = nps_df,
  ps = ps_df,
  acs_pop = acs_pop,
  mod = mod_int,
  adjust = TRUE,
  bootstrap = TRUE,
  L = 100
)
elapsed <- Sys.time() - start_time
cat("✓ Completed in", round(elapsed, 2), attr(elapsed, "units"), "\n")

mrp_int_r_boot <- mrp_int_boot$puma_summary_mrpr_bootstrap

# ============================================================
# COMBINE AND ANALYZE RESULTS
# ============================================================
cat("\n=== Analyzing Results ===\n\n")

all_results <- bind_rows(
  direst,
  mrp_p_baseline,
  mrp_r_baseline,
  mrp_r_boot,
  mrp_int_p_baseline,
  mrp_int_r_baseline,
  mrp_int_r_boot
)

# Calculate performance metrics
int_score <- function(alpha, y, lower, upper) {
  (upper - lower) +
    (2 / alpha) * (lower - y) * (y < lower) +
    (2 / alpha) * (y - upper) * (y > upper)
}

performance_metrics <- true_values %>%
  left_join(all_results, by = "PUMA") %>%
  group_by(model) %>%
  summarize(
    MSE = mean((response_true - point_est)^2),
    MAB = mean(abs(response_true - point_est)),
    Coverage = mean(between(response_true, lower_CI, upper_CI)),
    `Int. Score` = mean(int_score(alpha, response_true, lower_CI, upper_CI)),
    `Avg Width` = mean(upper_CI - lower_CI),
    .groups = "drop"
  ) %>%
  arrange(match(model, c("direst", "mrp-p", "mrp-r", "mrp-r-bootstrap",
                         "mrp-p-int", "mrp-r-int", "mrp-r-int-bootstrap")))

cat("=== Overall Performance Metrics ===\n")
print(performance_metrics, width = 120)

# Compare baseline vs bootstrap
cat("\n=== Bootstrap Impact ===\n\n")

cat("MRP-R:\n")
mrp_compare <- performance_metrics %>%
  filter(model %in% c("mrp-r", "mrp-r-bootstrap")) %>%
  select(model, Coverage, MSE, MAB, `Int. Score`, `Avg Width`)
print(mrp_compare)

coverage_improvement_mrp <- mrp_compare %>%
  filter(model == "mrp-r-bootstrap") %>%
  pull(Coverage) -
  mrp_compare %>%
  filter(model == "mrp-r") %>%
  pull(Coverage)

cat("\n  → Coverage improvement:", sprintf("+%.1f%%", coverage_improvement_mrp * 100), "\n")

cat("\nMRP-INT-R:\n")
mrp_int_compare <- performance_metrics %>%
  filter(model %in% c("mrp-r-int", "mrp-r-int-bootstrap")) %>%
  select(model, Coverage, MSE, MAB, `Int. Score`, `Avg Width`)
print(mrp_int_compare)

coverage_improvement_int <- mrp_int_compare %>%
  filter(model == "mrp-r-int-bootstrap") %>%
  pull(Coverage) -
  mrp_int_compare %>%
  filter(model == "mrp-r-int") %>%
  pull(Coverage)

cat("\n  → Coverage improvement:", sprintf("+%.1f%%", coverage_improvement_int * 100), "\n")

cat("\n✅ All tests complete\n")
