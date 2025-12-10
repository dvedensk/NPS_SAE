library(readr)
library(dplyr)
library(sampling)
source(file.path("code", "sampling_functions.R"))

set.seed(99)

# Load population file
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv"))
# Note: Keep AGEP numeric for sampling functions that use it in weights
# Convert categorical variables to factors for modeling

# Draw test samples using finalized default weights
# PS: WAGP=0.05, PWGTP=-0.2, samp_frac=0.005 → DDC ≈ 0.0007 (essentially unbiased) for PUBCOV
ps1 <- get_strat_PS(
  pop_df = acs_pop,
  samp_frac = .005,
  weight_config = list(WAGP = 0.05, PWGTP = -0.2)
)

# NPS: PWGTP=0.3, AGEP=0.7, internet_only=FALSE, samp_frac=0.05 → DDC ≈ -0.076 (strong bias) for PUBCOV
# Effective sample size: ESS ≈ 5.2 (with n ≈ 73,000)
nps <- get_NPS(
  pop_df = acs_pop,
  samp_frac = .05,
  weight_config = list(PWGTP = 0.3, AGEP = 0.7)
  # internet_only defaults to FALSE
)

# Extract data frames and convert categorical variables to factors for modeling
ps_df <- acs_pop[ps$idx, ] %>%
  mutate(
    AGEP_binned = factor(AGEP_binned), # Use binned version for modeling
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

nps_df <- acs_pop[nps$idx, ] %>%
  mutate(
    AGEP_binned = factor(AGEP_binned), # Use binned version for modeling
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# Save test samples to /data
save(ps, nps, ps_df, nps_df, file = "data/test_sample.RData")
cat("\n-------------------------------------------------\n")
cat("\n Test samples saved to data/test_sample.RData!\n")
