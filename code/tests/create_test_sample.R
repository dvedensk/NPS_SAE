library(readr)
library(dplyr)
library(sampling)
source(file.path("code", "sampling_functions.R"))

set.seed(99)

# Load population file
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv")) %>%
  mutate(
    AGEP = factor(AGEP),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# Draw test samples
# PS: Optimized weights for near-zero correlation with HICOV
ps <- get_strat_PS(pop_df = acs_pop, samp_frac = .002)

# NPS: Optimized weights for small positive correlation with HICOV
# Uses default optimized weights (w_cit=5, w_wage=-10, w_pwgt=-25)
# and reduced sample fraction (5% vs 20%) for stronger selection bias
nps <- get_NPS(
  pop_df = acs_pop,
  samp_frac = .05,
  include_internet = FALSE
)

# Extract data frames
ps_df <- acs_pop[ps$idx, ]
nps_df <- acs_pop[nps$idx, ]

# Save test samples to /data
save(ps, nps, ps_df, nps_df, file = "data/test_sample.RData")
cat("Test samples saved to data/test_sample.RData\n")
