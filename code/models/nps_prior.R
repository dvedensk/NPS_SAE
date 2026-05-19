# Helper function; take draws from Stan, 
# get draws from the posterior predictive, 
# and aggregate into area-level predictions
format_stan_output <- function(
  X,
  beta_draws, 
  eta_draws,
  domain_index,
  PUMA, 
  PUMA_levels, 
  typeIerr, 
  model_name
) {
  n <- nrow(X)
  S <- nrow(beta_draws)

  # Posterior predictive for PS units
  eta <- tcrossprod(X, beta_draws) # n x S
  if (!is.null(eta_draws)) {
    eta_domain <- t(eta_draws[, domain_index, drop = FALSE])
    eta <- eta + eta_domain
  }
  P <- plogis(eta)
  U <- matrix(runif(n * S), n, S)
  y_pred <- (U < P)
  storage.mode(y_pred) <- "integer"

  # Area-level summaries (robust to factor/character PUMA)
  group_idx <- as.integer(factor(PUMA, levels = PUMA_levels))
  by_PUMA_sums <- matrix(0L, nrow = length(PUMA_levels), ncol = S)
  for (i in seq_len(length(group_idx))) {
    g <- group_idx[i]
    if (!is.na(g)) by_PUMA_sums[g, ] <- by_PUMA_sums[g, ] + y_pred[i, ]
  }
  denom <- tabulate(group_idx, nbins = length(PUMA_levels))
  denom[denom == 0] <- NA_integer_
  by_PUMA_means <- sweep(by_PUMA_sums, 1, denom, "/")

  point_est <- rowMeans(by_PUMA_means, na.rm = TRUE)
  CI <- matrixStats::rowQuantiles(
    by_PUMA_means,
    probs = c(typeIerr / 2, 1 - typeIerr / 2),
    na.rm = TRUE
  )

  res_df <- data.frame(
    PUMA = PUMA_levels,
    point_est = point_est,
    lower_CI = CI[, 1],
    upper_CI = CI[, 2],
    model = model_name,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  res_df
}

scale_weight_covariate <- function(x, scale_covariate) {
  # Keeps optional z-scoring for weight-as-covariate usage; helps match standardized priors
  x <- as.numeric(x)
  if (!scale_covariate) {
    return(x)
  }
  mu <- mean(x)
  sigma <- stats::sd(x)
  if (is.na(sigma) || sigma == 0) {
    return(x - mu)
  }
  (x - mu) / sigma
}

prepare_weight_covariate <- function(weight_covariate, n, n_np, scale_covariate) {
  if (is.null(weight_covariate)) {
    return(list(ps = NULL, nps = NULL))
  }
  if (!is.list(weight_covariate)) {
    stop("weight_covariate must be a list with components ps and nps.")
  }
  ps <- weight_covariate$ps
  nps <- weight_covariate$nps
  if (is.null(ps) || is.null(nps)) {
    stop("weight_covariate must include both ps and nps vectors.")
  }
  if (length(ps) != n || length(nps) != n_np) {
    stop("weight_covariate lengths must match y and y_NP lengths.")
  }

  list(
    ps = scale_weight_covariate(ps, scale_covariate),
    nps = scale_weight_covariate(nps, scale_covariate)
  )
}

calculate_adaptive_power_prior_a <- function(y_ps, X_ps, y_nps, X_nps, verbose = TRUE) {
  # Fit GLMs to get MLEs
  glm_ps <- glm(y_ps ~ X_ps - 1, family = binomial())
  glm_nps <- glm(y_nps ~ X_nps - 1, family = binomial())

  beta_ps_hat <- coef(glm_ps)
  beta_nps_hat <- coef(glm_nps)
  diff_vec <- beta_ps_hat - beta_nps_hat

  # Hotelling T^2 test
  n_ps <- length(y_ps)
  p_coef <- ncol(X_ps)
  cov_beta_ps <- vcov(glm_ps)

  # T^2 statistic
  t2 <- as.numeric(t(diff_vec) %*% solve(cov_beta_ps) %*% diff_vec)

  # Convert to F-statistic
  f_scaling <- p_coef * (n_ps - 1) / (n_ps - p_coef)
  f_stat <- t2 / f_scaling

  # p-value becomes the power prior exponent
  adaptive_a <- pf(f_stat, p_coef, n_ps - p_coef, lower.tail = FALSE)

  if (verbose) {
    cat(
      "\nAdaptive Power Prior Exponent Calculation:",
      "\n- Hotelling T^2:", round(t2, 4),
      "\n- F-statistic:", round(f_stat, 4),
      "\n- p-value:", round(adaptive_a, 4),
      "\n- Adaptive a:", round(adaptive_a, 4),
      "\n\nInterpretation:",
      "\n- a ≈ 1.0: Strong PS-NPS agreement → NPS gets high weight",
      "\n- a ≈ 0.5: Moderate agreement → NPS gets moderate weight",
      "\n- a ≈ 0.0: Poor agreement → NPS gets low weight (high bias suspected)",
      "\n\n"
    )
  }

  return(list(
    a = adaptive_a,
    t2 = t2,
    f_stat = f_stat,
    p_value = adaptive_a,
    beta_ps = beta_ps_hat,
    beta_nps = beta_nps_hat,
    diff = diff_vec
  ))
}

# Main function for getting predictions from all four prior specifications
nps_prior_mcmc <- function(
    d_mod, 
    # md_mod
    pp_mod,
    y, 
    X, 
    y_NP, 
    X_NP, 
    wts, 
    PUMA, 
    PUMA_levels, 
    raking_or_pl,
    typeIerr = 0.05, 
    niter = 4000, 
    warmup = 1000,
    chains = 4,
    seed = 99,
    which_prior = c("pp", "d", "dl", "dl10"),
    a = NULL,
    domain_ps = NULL,
    domain_nps = NULL,
    domain_levels = NULL,
    weight_covariate = NULL,
    scale_weight_covariate = TRUE,
    refresh = 10,
    threads_per_chain = NULL,
    parallel_chains = NULL
) {
  n <- length(y)
  n_np <- length(y_NP)

  which_prior <- unique(which_prior)
  valid_priors <- c("pp", "d", "dl", "dl10")
  if (!all(which_prior %in% valid_priors)) {
    stop("which_prior must be one or more of: pp, d, dl, dl10.")
  }

  if(is.null(wts)) { # If not using pseudolikelihood
    wts <- rep(1, n)
  } else { # If using pseudolikelihood
    stopifnot(
        "must provide one weight for every observation" = 
            length(wts) == n
    )
    # Rescale the weights to sum to the sample size
    if(sum(wts) != n) {
        wts <- n * wts / sum(wts)
    }
  }

  X_eff <- X
  X_NP_eff <- X_NP

  weight_covs <- prepare_weight_covariate(weight_covariate, n, n_np, scale_weight_covariate)
  if (!is.null(weight_covs$ps)) {
    X_eff <- cbind(X_eff, weight_covs$ps)
    X_NP_eff <- cbind(X_NP_eff, weight_covs$nps)
  }

  p <- ncol(X_eff)

  # TODO: regression test after renaming md -> d
  # TODO: implement `mixed-distance`
  need_d <- any(which_prior %in% c("d", "dl", "dl10"))
  need_pp <- "pp" %in% which_prior
  need_glm <- need_d || (need_pp && is.null(a))

  if (need_glm) {
    # Get MLEs for regression coefficients from each sample
    # PS
    glm_fit <- glm(y ~ X_eff - 1, binomial())
    beta_hat <- coef(glm_fit) %>% 
        as.numeric()

    stopifnot("Some coefficients are NA; X may be ill-conditioned." = !anyNA(beta_hat))

    # NPS
    glm_fit_NP <- glm(y_NP ~ X_NP_eff - 1, binomial())
    beta_NP_hat <- coef(glm_fit_NP) %>% 
        as.numeric()

    stopifnot("Some coefficients are NA; X_NP may be ill-conditioned." = !anyNA(beta_NP_hat))
    
    diff_vec <- beta_hat - beta_NP_hat
    sq_dist <- diff_vec^2
  }

  # Set up basic Stan data list
  if (is.null(domain_ps) || is.null(domain_nps)) {
    stop("domain_ps and domain_nps must be provided for domain random effects.")
  }
  if (is.null(domain_levels)) {
    domain_levels <- sort(unique(c(domain_ps, domain_nps)))
  }
  J <- length(domain_levels)

  stan_data <- list(
    n = n,
    p = p,
    J = J,
    X = X_eff,
    y = as.integer(y),
    w = as.numeric(wts),
    domain = as.integer(domain_ps)
  )

  run_sample <- function(mod, data) {
    sample_args <- list(
      data = data,
      iter_sampling = niter,
      iter_warmup = warmup,
      chains = chains,
      seed = seed,
      refresh = refresh
    )
    if (!is.null(threads_per_chain)) {
      sample_args$threads_per_chain <- threads_per_chain
    }
    if (!is.null(parallel_chains)) {
      sample_args$parallel_chains <- parallel_chains
    }
    do.call(mod$sample, sample_args)
  }

  if (need_d) {
    d_base_data <- c(stan_data, list(
      beta_prior_mean = beta_NP_hat
    ))
  }

  if ("d" %in% which_prior) {
    # Distance prior
    prior_sds <- sqrt(sq_dist)

    d_data <- c(d_base_data, list(
      beta_prior_sd = prior_sds
    ))

    d_fit <- run_sample(d_mod, d_data)

    print(d_fit$summary()[,"rhat"])

    # Extract posterior draws
    beta_draws <- d_fit$draws("beta", format = "matrix")
    eta_draws <- d_fit$draws("eta_domain", format = "matrix")
    # Format into area-level summaries
    d_res <- format_stan_output(
      X_eff, 
      beta_draws, 
      eta_draws,
      domain_ps,
      PUMA, 
      PUMA_levels, 
      typeIerr, 
      paste0("NPS Prior w/ ", raking_or_pl, ": Distance")
    )
  }

  if ("dl" %in% which_prior || "dl10" %in% which_prior) {
    # for distance-log priors, we estimate the variance of the MLE 
    # of each coefficient based on the NPS 
    NP_vars <- glm_fit_NP %>% 
      vcov() %>% 
      diag() %>% 
      as.numeric()

    maxes <- pmax(sq_dist, NP_vars) 
  }

  if ("dl" %in% which_prior) {
    # Distance-log
    prior_sds <- sqrt(maxes / log(n_np))

    dl_data <- c(d_base_data, list(
      beta_prior_sd = prior_sds
    ))

    dl_fit <- run_sample(d_mod, dl_data)

    print(dl_fit$summary()[,"rhat"])

    # Extract posterior draws
    beta_draws <- dl_fit$draws("beta", format = "matrix")
    eta_draws <- dl_fit$draws("eta_domain", format = "matrix")
    # Format into area-level summaries
    dl_res <- format_stan_output(
      X_eff, 
      beta_draws, 
      eta_draws,
      domain_ps,
      PUMA, 
      PUMA_levels, 
      typeIerr, 
      paste0("NPS Prior w/ ", raking_or_pl, ": Distance-log")
    )
  }

  if ("dl10" %in% which_prior) {
    prior_sds <- sqrt(maxes / log10(n_np))

    dl_ten_data <- c(d_base_data, list(
      beta_prior_sd = prior_sds
    ))

    dl_ten_fit <- run_sample(d_mod, dl_ten_data)

    print(dl_ten_fit$summary()[,"rhat"])

    # Extract posterior draws
    beta_draws <- dl_ten_fit$draws("beta", format = "matrix")
    eta_draws <- dl_ten_fit$draws("eta_domain", format = "matrix")
    # Format into area-level summaries
    dl_ten_res <- format_stan_output(
      X_eff, 
      beta_draws, 
      eta_draws,
      domain_ps,
      PUMA, 
      PUMA_levels, 
      typeIerr, 
      paste0("NPS Prior w/ ", raking_or_pl, ": Distance-log10")
    )
  }

  if (need_pp && is.null(a)) {
    # Power Prior
    # Set power prior exponent to p-value corresponding to Hotelling T^2 test
    cov_beta_hat <- vcov(glm_fit)
    t2 <- crossprod(diff_vec, (cov_beta_hat %>% chol() %>% chol2inv()) %*% diff_vec)
    t2 <- as.numeric(t2)
    f = p * (n - 1) / (n - p)
    Fstat = t2 / f
    a <- pf(Fstat, p, n - p, lower.tail = F)
  }

  if (need_pp) {
    pp_data <- c(stan_data, list(
      n_np = n_np,
      X_np = X_NP_eff, 
      y_np = as.integer(y_NP), 
      a = a,
      domain_np = as.integer(domain_nps)
    ))

    pp_fit <- run_sample(pp_mod, pp_data)

    print(pp_fit$summary()[,"rhat"])

    # Extract posterior draws
    beta_draws <- pp_fit$draws("beta", format = "matrix")
    eta_draws <- pp_fit$draws("eta_domain", format = "matrix")
    # Format into area-level summaries
    pp_res <- format_stan_output(
      X_eff, 
      beta_draws, 
      eta_draws,
      domain_ps,
      PUMA, 
      PUMA_levels, 
      typeIerr, 
      paste0("NPS Prior w/ ", raking_or_pl, ": Power Prior")
    )
  }

  res_list <- list()
  if ("d" %in% which_prior) res_list[["d"]] <- d_res
  if ("dl" %in% which_prior) res_list[["dl"]] <- dl_res
  if ("dl10" %in% which_prior) res_list[["dl10"]] <- dl_ten_res
  if ("pp" %in% which_prior) res_list[["pp"]] <- pp_res
  res_df <- do.call(rbind, res_list)
  
  return(res_df)
}
