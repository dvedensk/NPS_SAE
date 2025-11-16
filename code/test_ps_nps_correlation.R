library(dplyr)
library(readr)
library(sampling)

source("~/dev/MSS_NPS/code/sampling_functions.R")

# Load population data
acs_pop <- read_csv("~/dev/MSS_NPS/data/ACS_NPS_pop.csv")

# Fixed weights
PS_WEIGHTS <- list(w_cit = 8, w_wage = 8, w_pwgt = -6)
NPS_WEIGHTS <- list(w_cit = 25, w_wage = -10, w_pwgt = -25)

# Sample fractions
PS_SAMP_FRAC <- 0.002
NPS_SAMP_FRAC <- 0.05

# Number of replications
N_REPS <- 30

cat("Testing correlations with HICOV:\n")
cat("PS weights: w_cit =", PS_WEIGHTS$w_cit, ", w_wage =", PS_WEIGHTS$w_wage, ", w_pwgt =", PS_WEIGHTS$w_pwgt, "\n")
cat("NPS weights: w_cit =", NPS_WEIGHTS$w_cit, ", w_wage =", NPS_WEIGHTS$w_wage, ", w_pwgt =", NPS_WEIGHTS$w_pwgt, "\n")
cat("PS sample fraction:", PS_SAMP_FRAC, "\n")
cat("NPS sample fraction:", NPS_SAMP_FRAC, "\n")
cat("Running", N_REPS, "replications...\n\n")

ps_correlations <- numeric(N_REPS)
nps_correlations <- numeric(N_REPS)

for (i in 1:N_REPS) {
  set.seed(i)

  # Draw PS sample
  ps <- get_strat_PS(
    pop_df = acs_pop,
    samp_frac = PS_SAMP_FRAC,
    w_cit = PS_WEIGHTS$w_cit,
    w_wage = PS_WEIGHTS$w_wage,
    w_pwgt = PS_WEIGHTS$w_pwgt
  )

  # Draw NPS sample
  nps <- get_NPS(
    pop_df = acs_pop,
    samp_frac = NPS_SAMP_FRAC,
    w_cit = NPS_WEIGHTS$w_cit,
    w_wage = NPS_WEIGHTS$w_wage,
    w_pwgt = NPS_WEIGHTS$w_pwgt,
    include_internet = FALSE
  )

  # Create indicator vectors
  in_ps <- rep(0, nrow(acs_pop))
  in_nps <- rep(0, nrow(acs_pop))
  in_ps[ps$idx] <- 1
  in_nps[nps$idx] <- 1

  # Calculate correlations with HICOV
  ps_correlations[i] <- cor(in_ps, acs_pop$HICOV)
  nps_correlations[i] <- cor(in_nps, acs_pop$HICOV)

  if (i %% 10 == 0) {
    cat("Completed", i, "replications...\n")
  }
}

cat("\n=== PS Sample Correlation with HICOV ===\n")
cat("Mean correlation:", mean(ps_correlations), "\n")
cat("Median correlation:", median(ps_correlations), "\n")
cat("SD correlation:", sd(ps_correlations), "\n")
cat("Min correlation:", min(ps_correlations), "\n")
cat("Max correlation:", max(ps_correlations), "\n")

cat("\n=== NPS Sample Correlation with HICOV ===\n")
cat("Mean correlation:", mean(nps_correlations), "\n")
cat("Median correlation:", median(nps_correlations), "\n")
cat("SD correlation:", sd(nps_correlations), "\n")
cat("Min correlation:", min(nps_correlations), "\n")
cat("Max correlation:", max(nps_correlations), "\n")
