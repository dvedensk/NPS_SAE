# Stratified Probability Sample (PS) by PUMA
# Recommended samp_frac: 0.005 (0.5% → ~7,378 individuals, 10x smaller than NPS)
# Default weights: w_wage=0.5, w_pwgt=-0.5
# These achieve DDC ≈ 0.0028 with PUBCOV (essentially unbiased)
# Coverage: 99.9% of 281 PUMAs (280.9 on average)
# Zero-variance PUMAs: ~4 (1.4%) - acceptable for small area estimation
get_strat_PS <- function(pop_df, samp_frac, w_wage = 0.5, w_pwgt = -0.5, w_povpip = 0, w_ssp = 0, w_agep = 0) {
  PUMAs <- unique(pop_df$PUMA)
  weights <- idx <- c()
  for (PUMA in PUMAs) {
    puma_ids <- which(pop_df$PUMA == PUMA)
    ps <- get_PS(
      pop_df = pop_df[puma_ids, ],
      samp_frac = samp_frac,
      w_wage = w_wage,
      w_pwgt = w_pwgt,
      w_povpip = w_povpip,
      w_ssp = w_ssp,
      w_agep = w_agep
    )
    weights <- c(weights, ps$weights)
    idx <- c(idx, ps$idx)
  }
  # Note: This check was previously implemented incorrectly (only checked the last PUMA in the loop).
  if (sum(weights == 1) > 0) {
    stop("Some weights were equal to 1. Adjust the informative sampling weights.")
  }
  return(list(weights = weights, idx = idx))
}

get_PS <- function(pop_df,
                   # type = "PPS", # seems like this input is no longer used
                   samp_frac = .002,
                   w_wage = 0.5,
                   w_pwgt = -0.5,
                   w_povpip = 0,
                   w_ssp = 0,
                   w_agep = 0) { # type is PPS or some other size variable
  # sets minimal sample size to 1 (higher minimums resulted in the `sum(weights == 1) > 0` check to go off
  sample_size <- max(floor(nrow(pop_df) * samp_frac), 1)
  # create the size variable that will be sampled in proportion to
  # (making this variable a function of the survey weights induces an informative design)
  # All variables are scaled for comparability
  size_var <- as.numeric(exp(w_wage * scale(pop_df$WAGP) +
    w_pwgt * scale(pop_df$PWGTP) +
    w_povpip * scale(pop_df$POVPIP) +
    w_ssp * scale(pop_df$SSP) +
    w_agep * scale(pop_df$AGEP)))
  inclusion_probs <- inclusionprobabilities(size_var, sample_size)
  inclusion_probs <- inclusion_probs / sum(inclusion_probs) * sample_size
  # survey weights are inverse probabilities of selection
  weights <- 1 / inclusion_probs
  # Draw Poisson sample: https://en.wikipedia.org/wiki/Poisson_sampling
  sample_idx <- which(UPpoisson(inclusion_probs) == 1)
  sample_size <- length(sample_idx)
  weights <- weights[sample_idx]

  return(list(weights = weights, idx = sample_idx))
}

# Non-Probability Sample (NPS)
# Recommended samp_frac: 0.05 (5% → ~73,781 individuals, 10x larger than PS)
# Default weights: w_pwgt=0.2, w_agep=0.8, internet_only=FALSE
# These achieve DDC ≈ -0.093 with PUBCOV (strong selection bias for testing bias correction)
# Coverage: 100% of 281 PUMAs (perfect domain coverage)
# Effective sample size: n_eff ≈ 6 (with actual n ≈ 74,000)
# Alternative configs:
#   - For stronger bias: w_agep=1.0 → DDC=-0.123, but n_eff=4 (extreme)
#   - For moderate bias: w_pwgt=0.3, w_agep=0.7 → DDC=-0.048, n_eff=22
#   - With internet_only=TRUE: DDC ≈ -0.0014 (realistic but weak bias)
# See NPS_WEIGHT_SENSITIVITY.md for full trade-off analysis
get_NPS <- function(pop_df,
                    w_wage = 0,
                    w_pwgt = 0.2,
                    w_povpip = 0,
                    w_ssp = 0,
                    w_agep = 0.8,
                    samp_frac = .05,
                    internet_only = FALSE) {
  if (internet_only) {
    pop_df <- filter(pop_df, ACCESSINET < 3) # 3 indicates no access
  }

  # All variables are scaled for comparability
  size_var <- as.numeric(exp(w_wage * scale(pop_df$WAGP) +
    w_pwgt * scale(pop_df$PWGTP) +
    w_povpip * scale(pop_df$POVPIP) +
    w_ssp * scale(pop_df$SSP) +
    w_agep * scale(pop_df$AGEP)))

  # Simple inclusion probability calculation (preserves extreme bias)
  # This approach achieves much higher DDC than inclusionprobabilities()
  inclusion_probs <- samp_frac * size_var / mean(size_var)
  inclusion_probs <- pmin(inclusion_probs, 1) # Cap at 1.0
  inclusion_probs[is.na(inclusion_probs)] <- 0

  weights <- 1 / inclusion_probs
  weights[inclusion_probs == 0] <- 0

  sample_idx <- which(UPpoisson(inclusion_probs) == 1)
  sample_size <- length(sample_idx)
  weights <- weights[sample_idx]

  # for assessing data defect index, see test_ps_nps_correlation.R

  # return weights since we may use them for calculating properties of the sample, but
  # we will not use them as...
  return(list(weights = weights, idx = sample_idx))
}
