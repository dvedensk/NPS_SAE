library(readr)
library(dplyr)
library(sampling)
source(file.path("code", "sampling_functions.R"))

set.seed(99)

# Load population file
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv"))
# Note: Keep AGEP numeric for sampling functions that use it in weights
# Convert to factor only after sampling for modeling

# Draw test samples using finalized default weights
# PS: w_wage=0.5, w_pwgt=-0.5, samp_frac=0.005 → DDC ≈ 0.0028 (essentially unbiased)
ps <- get_strat_PS(pop_df = acs_pop, samp_frac = .005)

# NPS: w_pwgt=0.2, w_agep=0.8, internet_only=FALSE, samp_frac=0.05 → DDC ≈ -0.093 (strong bias)
# Effective sample size: ESS ≈ 3.2 (with n ≈ 74,000)
nps <- get_NPS(
  pop_df = acs_pop,
  samp_frac = .05
  # Uses defaults: w_pwgt=0.2, w_agep=0.8, internet_only=FALSE
)

# Extract data frames and convert categorical variables to factors for modeling
ps_df <- acs_pop[ps$idx, ] %>%
  mutate(
    AGEP = factor(AGEP_binned),  # Use binned version for modeling
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

nps_df <- acs_pop[nps$idx, ] %>%
  mutate(
    AGEP = factor(AGEP_binned),  # Use binned version for modeling
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# Save test samples to /data
save(ps, nps, ps_df, nps_df, file = "data/test_sample.RData")
cat("Test samples saved to data/test_sample.RData\n")
