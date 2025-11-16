get_strat_PS <- function(pop_df, samp_frac, w_wage = 0.5, w_pwgt = -0.5, w_cit = 0.5) {
  PUMAs <- unique(pop_df$PUMA)
  weights <- idx <- c()
  for (PUMA in PUMAs) {
    puma_ids <- which(pop_df$PUMA == PUMA)
    ps <- get_PS(
      pop_df = pop_df[puma_ids, ],
      samp_frac = samp_frac,
      w_wage = w_wage,
      w_pwgt = w_pwgt,
      w_cit = w_cit
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
                   w_cit = 0.5,
                   w_wage = 0.5,
                   w_pwgt = -0.5) { # type is PPS or some other size variable
  # sets minimal sample size to 1 (higher minimums resulted in the `sum(weights == 1) > 0` check to go off
  sample_size <- max(floor(nrow(pop_df) * samp_frac), 1)
  # create the size variable that will be sampled in proportion to
  # (making this variable a function of the survey weights induces an informative design)
  # All variables are scaled for comparability
  size_var <- as.numeric(exp(w_wage * scale(pop_df$WAGP) +
    w_pwgt * scale(pop_df$PWGTP) +
    w_cit * scale(pop_df$CIT == 5)))
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

get_NPS <- function(pop_df,
                    w_cit = 25,
                    w_wage = -10,
                    w_pwgt = -25,
                    samp_frac = .05,
                    include_internet = F) {
  if (!include_internet) {
    pop_df <- filter(pop_df, ACCESSINET < 3) # 3 indicates no access
  }
  sample_size <- floor(nrow(pop_df) * samp_frac)
  # if R \propto G, ddi = 1, so start there and add noise
  # All variables are scaled for comparability
  size_var <- as.numeric(exp(w_wage * scale(pop_df$WAGP) +
    w_pwgt * scale(pop_df$PWGTP) +
    w_cit * scale(pop_df$CIT == 5)))
  inclusion_probs <- inclusionprobabilities(size_var, sample_size)
  inclusion_probs <- inclusion_probs / sum(inclusion_probs) * sample_size
  weights <- 1 / inclusion_probs
  sample_idx <- which(UPpoisson(inclusion_probs) == 1)
  sample_size <- length(sample_idx)
  weights <- weights[sample_idx]

  # for assessing data defect index, see test_ps_nps_correlation.R

  # return weights since we may use them for calculating properties of the sample, but
  # we will not use them as...
  return(list(weights = weights, idx = sample_idx))
}
