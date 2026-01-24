# Stratified Probability Sample (PS) by PUMA
# Recommended samp_frac: 0.005 (0.5% → ~7,264 individuals, 10x smaller than NPS)
# Default weights: WAGP=0.05, PWGTP=-0.2 (minimal informative sampling)
# These achieve DDC ≈ 0.0007 with PUBCOV (essentially unbiased)
# ESS: ~1,367,527 (median: 15,291) - extremely high effective sample size
# Coverage: 100% of 281 PUMAs (perfect stratification)
# Sample bias: +0.46 pp (60.8% vs pop 60.3%) - minimal bias
# Zero-variance PUMAs: ~0.1 (0.0%) - essentially zero
#
# weight_config: named list specifying weight variables and coefficients
#   e.g., list(PWGTP = 0.1, AGEP = 0.5) or list(WAGP = 0.05, PWGTP = -0.2)
get_strat_PS <- function(pop_df, samp_frac, weight_config = list(WAGP = 0.05, PWGTP = -0.2)) {
  PUMAs <- unique(pop_df$PUMA)
  weights <- idx <- c()
  for (PUMA in PUMAs) {
    puma_ids <- which(pop_df$PUMA == PUMA)
    ps <- get_PS(
      pop_df = pop_df[puma_ids, ],
      samp_frac = samp_frac,
      weight_config = weight_config
    )
    weights <- c(weights, ps$weights)
    # Map subset indices back to original dataframe indices
    idx <- c(idx, puma_ids[ps$idx])
  }
  # Note: This check was previously implemented incorrectly (only checked the last PUMA in the loop).
  if (sum(weights == 1) > 0) {
    stop("Some weights were equal to 1. Adjust the informative sampling weights.")
  }
  return(list(weights = weights, idx = idx))
}

get_PS <- function(pop_df,
                   samp_frac = .002,
                   weight_config = list(PWGTP = 0.1)) {
  # sets minimal sample size to 1 (higher minimums resulted in the `sum(weights == 1) > 0` check to go off
  sample_size <- max(floor(nrow(pop_df) * samp_frac), 1)

  # create the size variable that will be sampled in proportion to
  # (making this variable a function of the survey weights induces an informative design)
  # All variables are scaled for comparability
  # Build size_var dynamically from weight_config
  linear_terms <- 0
  for (var_name in names(weight_config)) {
    if (!var_name %in% names(pop_df)) {
      stop(paste0("Variable '", var_name, "' not found in pop_df"))
    }
    linear_terms <- linear_terms + weight_config[[var_name]] * scale(pop_df[[var_name]])
  }
  size_var <- as.numeric(exp(linear_terms))

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
# Recommended samp_frac: 0.05 (5% → ~73,000 individuals, 10x larger than PS)
# Default weights: PWGTP=0.3, AGEP=0.7 (for PUBCOV response)
# These achieve DDC ≈ -0.076 with PUBCOV (strong selection bias for testing bias correction)
# Coverage: 100% of 281 PUMAs (perfect domain coverage)
# Effective sample size: ESS ≈ 5.2 (with actual n ≈ 73,000)
# Alternative configs:
#   - For stronger bias: PWGTP=0.2, AGEP=0.8 → DDC≈-0.093, ESS≈6
#   - For extreme bias: PWGTP=0, AGEP=1.0 → DDC≈-0.123, ESS≈4
#   - With internet_only=TRUE: DDC ≈ -0.0014 (realistic but weak bias)
# See NPS_WEIGHT_SENSITIVITY.md for full trade-off analysis
#
# weight_config: named list specifying weight variables and coefficients
#   e.g., list(PWGTP = 0.3, AGEP = 0.7) or list(WAGP = 0.5)
get_NPS <- function(pop_df,
                    samp_frac = .05,
                    weight_config = list(PWGTP = 0.3, AGEP = 0.7),
                    internet_only = FALSE,
                    rescale_after_cap = FALSE) {
  if (internet_only) {
    pop_df <- filter(pop_df, ACCESSINET < 3) # 3 indicates no access
  }

  # All variables are scaled for comparability
  # Build size_var dynamically from weight_config
  linear_terms <- 0
  for (var_name in names(weight_config)) {
    if (!var_name %in% names(pop_df)) {
      stop(paste0("Variable '", var_name, "' not found in pop_df"))
    }
    linear_terms <- linear_terms + weight_config[[var_name]] * scale(pop_df[[var_name]])
  }
  size_var <- as.numeric(exp(linear_terms))

  # Simple inclusion probability calculation (preserves extreme bias)
  # This approach achieves much higher DDC than inclusionprobabilities()
  if (rescale_after_cap) {
    target <- samp_frac
    f_mean_prob <- function(k) {
      mean(pmin(k * size_var, 1), na.rm = TRUE) - target
    }
    upper <- samp_frac / mean(size_var, na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) {
      upper <- 1
    }
    iter <- 0
    while (f_mean_prob(upper) < 0 && iter < 50) {
      upper <- upper * 2
      iter <- iter + 1
    }
    if (f_mean_prob(upper) < 0) {
      inclusion_probs <- pmin(upper * size_var, 1)
    } else {
      k_star <- uniroot(f_mean_prob, lower = 0, upper = upper)$root
      inclusion_probs <- pmin(k_star * size_var, 1)
    }
  } else {
    inclusion_probs <- samp_frac * size_var / mean(size_var, na.rm = TRUE)
    inclusion_probs <- pmin(inclusion_probs, 1) # Cap at 1.0
  }
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

# Calculate Sample Diagnostics (DDC and ESS)
# Computes Data Defect Correlation (DDC) and Effective Sample Size (ESS) for both PS and NPS
#
# Arguments:
#   pop_df: Full population data frame
#   ps_sample: PS sample object from get_strat_PS() (list with idx and weights)
#   nps_sample: NPS sample object from get_NPS() (list with idx and weights)
#   response_var: Name of response variable (e.g., "PUBCOV", "HICOV")
#   print_results: If TRUE, prints diagnostic summary (default: FALSE)
#   sim_num: Simulation number for printing (optional)
#
# Returns:
#   List with ps_ddc, ps_ess, nps_ddc, nps_ess
#
# DDC measures selection bias (correlation between inclusion and response)
# ESS uses Meng's formula: ESS = (n_W/N) / (1 - n_W/N) * 1/rho^2
#   where n_W = n / (1 + CV_W^2) accounts for weight variability
calculate_sample_diagnostics <- function(pop_df,
                                         ps_sample,
                                         nps_sample,
                                         response_var,
                                         print_results = FALSE,
                                         sim_num = NULL) {
  # Create inclusion indicator vectors
  in_ps <- rep(0, nrow(pop_df))
  in_nps <- rep(0, nrow(pop_df))
  in_ps[ps_sample$idx] <- 1
  in_nps[nps_sample$idx] <- 1

  # Calculate DDC: correlation between inclusion indicator and response variable
  ps_ddc <- cor(in_ps, pop_df[[response_var]])
  nps_ddc <- cor(in_nps, pop_df[[response_var]])

  # Calculate ESS using Meng's formula: ESS = (n_W/N) / (1 - n_W/N) * 1/rho^2
  # Step 1: Calculate n_W = n / (1 + CV_W^2)

  # For PS (stratified): calculate per PUMA and sum
  ps_sample_temp <- pop_df[ps_sample$idx, ]
  ps_sample_temp$weight <- ps_sample$weights
  ps_n_W <- ps_sample_temp %>%
    group_by(PUMA) %>%
    summarise(
      n_puma = n(),
      CV_W = sd(weight) / mean(weight),
      n_W = n_puma / (1 + CV_W^2),
      .groups = "drop"
    ) %>%
    summarise(total_n_W = sum(n_W, na.rm = TRUE)) %>%
    pull(total_n_W)

  # For NPS (not stratified): calculate overall
  nps_sample_temp <- pop_df[nps_sample$idx, ]
  nps_sample_temp$weight <- nps_sample$weights
  nps_CV_W <- sd(nps_sample_temp$weight) / mean(nps_sample_temp$weight)
  nps_n_W <- length(nps_sample$idx) / (1 + nps_CV_W^2)

  # Step 2: ESS = (n_W/N) / (1 - n_W/N) * 1/rho^2
  ps_sf_adj <- ps_n_W / nrow(pop_df)
  ps_ess <- (ps_sf_adj / (1 - ps_sf_adj)) / (ps_ddc^2)

  nps_sf_adj <- nps_n_W / nrow(pop_df)
  nps_ess <- (nps_sf_adj / (1 - nps_sf_adj)) / (nps_ddc^2)

  # Print diagnostics if requested
  if (print_results) {
    if (!is.null(sim_num)) {
      cat("\n--- Sample Diagnostics (Simulation", sim_num, ") ---\n")
    } else {
      cat("\n--- Sample Diagnostics ---\n")
    }
    cat("PS:  n =", length(ps_sample$idx), "| DDC =", sprintf("%.4f", ps_ddc), "| ESS =", round(ps_ess), "\n")
    cat("NPS: n =", length(nps_sample$idx), "| DDC =", sprintf("%.4f", nps_ddc), "| ESS =", round(nps_ess), "\n\n")
  }

  return(list(
    ps_ddc = ps_ddc,
    ps_ess = ps_ess,
    nps_ddc = nps_ddc,
    nps_ess = nps_ess
  ))
}
