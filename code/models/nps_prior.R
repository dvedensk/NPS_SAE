# Helper function; take draws from Stan, 
# get draws from the posterior predictive, 
# and aggregate into area-level predictions
format_stan_output <- function(
  X,
  beta_draws, 
  PUMA, 
  PUMA_levels, 
  typeIerr, 
  model_name
) {
  n <- nrow(X)
  S <- nrow(beta_draws)

  # Posterior predictive for PS units
  # Draw linear predictor and then sample 1 with probability expit(eta), 0 otherwise
  # (Using latent Uniform(0,1) variates makes this extremely fast).
  eta <- tcrossprod(X, beta_draws) # n x S
  P <- plogis(eta)
  U <- matrix(runif(n * S), n, S)
  y_pred <- (U < P)
  storage.mode(y_pred) <- "integer"

  # Area-level summaries
  by_PUMA_sums <- rowsum(y_pred, group = PUMA, reorder = FALSE)
  by_PUMA_means <- by_PUMA_sums / as.numeric(table(PUMA))[PUMA_levels]
  point_est <- rowMeans(by_PUMA_means)
  CI <- rowQuantiles(
    by_PUMA_means,
    probs = c(typeIerr / 2, 1 - typeIerr / 2)
  )
  res_df <- data.frame(
    PUMA = PUMA_levels,
    point_est = point_est,
    lower_CI = CI[, 1],
    upper_CI = CI[, 2],
    model = model_name
  )

  return(res_df)
}

# Main function for getting predictions from all four prior specifications
nps_prior_mcmc <- function(
    md_mod, 
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
    seed = 99
) {
  n <- length(y)
  p <- ncol(X)
  n_NP <- length(y_NP)

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

  # Get MLEs for regression coefficients from each sample
  # PS
  glm_fit <- glm(y ~ X - 1, binomial())
  beta_hat <- coef(glm_fit) %>% 
      as.numeric()

  stopifnot("Some coefficients are NA; X may be ill-conditioned." = !anyNA(beta_hat))

  # NPS
  glm_fit_NP <- glm(y_NP ~ X_NP - 1, binomial())
  beta_NP_hat <- coef(glm_fit_NP) %>% 
      as.numeric()

  stopifnot("Some coefficients are NA; X_NP may be ill-conditioned." = !anyNA(beta_hat))
  
  diff_vec <- beta_hat - beta_NP_hat
  sq_dist <- diff_vec^2

  # Mixed-distance prior
  prior_sds <- sqrt(sq_dist)

  # Set up basic Stan data list
  stan_data <- list(
    n = n,
    p = p,
    X = X,
    y = y,
    w = as.numeric(wts)
  ) 

  md_base_data <- c(stan_data, list(
    beta_prior_mean = beta_NP_hat
  ))

  md_data <- c(md_base_data, list(
    beta_prior_sd = prior_sds
  ))

  md_fit <- md_mod$sample(
    data = md_data,
    iter_sampling = niter,
    iter_warmup = warmup,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 1,
    seed = seed,
    refresh = 10
  )

  rhat_md <- md_fit$summary()[,"rhat"]

  # Extract posterior draws
  beta_draws <- md_fit$draws("beta", format = "matrix")
  # Format into area-level summaries
  md_res <- format_stan_output(
    X, 
    beta_draws, 
    PUMA, 
    PUMA_levels, 
    typeIerr, 
    paste0("NPS Prior w/ ", raking_or_pl, ": Mixed-Distance")
  )

  # for mixed-distance-log priors, we estimate the variance of the MLE 
  # of each coefficient based on the NPS 
  NP_vars <- glm_fit_NP %>% 
    vcov() %>% 
    diag() %>% 
    as.numeric()

  maxes <- pmax(sq_dist, NP_vars) 

  # Mixed-distance-log
  prior_sds <- sqrt(maxes / log(n_NP))

  mdl_data <- c(md_base_data, list(
    beta_prior_sd = prior_sds
  ))

  mdl_fit <- md_mod$sample(
    data = mdl_data,
    iter_sampling = niter,
    iter_warmup = warmup,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 1,
    seed = seed,
    refresh = 10
  )

  rhat_mdl <- mdl_fit$summary()[,"rhat"]

  # Extract posterior draws
  beta_draws <- mdl_fit$draws("beta", format = "matrix")
  # Format into area-level summaries
  mdl_res <- format_stan_output(
    X, 
    beta_draws, 
    PUMA, 
    PUMA_levels, 
    typeIerr, 
    paste0("NPS Prior w/ ", raking_or_pl, ": Mixed-Distance-log")
  )

  prior_sds <- sqrt(maxes / log10(n_NP))

  mdl_ten_data <- c(md_base_data, list(
    beta_prior_sd = prior_sds
  ))

  mdl_ten_fit <- md_mod$sample(
    data = mdl_ten_data,
    iter_sampling = niter,
    iter_warmup = warmup,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 1,
    seed = seed,
    refresh = 10
  )

  rhat_mdl_ten <- mdl_ten_fit$summary()[,"rhat"]

  # Extract posterior draws
  beta_draws <- mdl_ten_fit$draws("beta", format = "matrix")
  # Format into area-level summaries
  mdl_ten_res <- format_stan_output(
    X, 
    beta_draws, 
    PUMA, 
    PUMA_levels, 
    typeIerr, 
    paste0("NPS Prior w/ ", raking_or_pl, ": Mixed-Distance-log10")
  )

  # Power Prior
  # Set power prior exponent to p-value corresponding to Hotelling T^2 test
  cov_beta_hat <- vcov(glm_fit)
  t2 <- crossprod(diff_vec, (cov_beta_hat %>% chol() %>% chol2inv()) %*% diff_vec)
  t2 <- as.numeric(t2)
  f = p * (n - 1) / (n - p)
  Fstat = t2 / f
  a <- pf(Fstat, p, n - p, lower.tail = F)

  pp_data <- c(stan_data, list(
    n_np = n_NP,
    X_np = X_NP, 
    y_np = y_NP, 
    a = a
  ))

  pp_fit <- pp_mod$sample(
    data = pp_data,
    iter_sampling = niter,
    iter_warmup = warmup,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 1,
    seed = seed,
    refresh = 10
  )

  rhat_pp <- pp_fit$summary()[,"rhat"]

  # Extract posterior draws
  beta_draws <- pp_fit$draws("beta", format = "matrix")
  # Format into area-level summaries
  pp_res <- format_stan_output(
    X, 
    beta_draws, 
    PUMA, 
    PUMA_levels, 
    typeIerr, 
    paste0("NPS Prior w/ ", raking_or_pl, ": Power Prior")
  )

  res_df <- rbind(
    md_res, 
    mdl_res, 
    mdl_ten_res,
    pp_res
  )

  save(rhat_md, rhat_mdl, rhat_mdl_ten, rhat_pp, file = paste0("rhat_", raking_or_pl, ".Rdata"))
  
  return(res_df)
}