#!/usr/bin/env Rscript
# Test MRP-R coverage with current implementation (weights, no bootstrap)

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

cat("\n=== Testing MRP-R Coverage (Baseline: weights, no bootstrap) ===\n\n")

# Compile Stan model
cat("Compiling Stan model...\n")
mod <- cmdstan_model(
  here("code/models", "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

# Run MRP with current settings
cat("\nRunning getMRP...\n")
start_time <- Sys.time()
mrp <- getMRP(
  MR = nps_df,
  ps = ps_df,
  acs_pop = acs_pop,
  bootstrap = FALSE,  # Baseline: no bootstrap
  L = 5
)
elapsed <- Sys.time() - start_time
cat("✅ Completed in", round(elapsed, 2), attr(elapsed, "units"), "\n\n")

# Extract results
mrpra <- mrp$puma_summary_mrpr
mrppb <- mrp$puma_summary_mrpp

# Calculate true values
response_var <- "PUBCOV"
alpha <- 0.05

true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    response_true = mean(.data[[response_var]], na.rm = TRUE),
    .groups = "drop"
  )

# Direct estimates
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

# Combine results
all_results <- bind_rows(direst, mrppb, mrpra)

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
    .groups = "drop"
  )

cat("\n=== Performance Metrics (Baseline) ===\n")
print(performance_metrics, width = 100)

cat("\n✅ Test complete\n")
