#' NPS Prior Method

#' Important: marginal posterior of beta is multivariate t with 
#' location `post_mn`, scale matrix `post_var * (n + 2 * a - 2) / (n + 2 * a)`, 
#' and degrees of freedom = `n + 2 * a`
#' Generate Monte Carlo samples from the marginal posterior 
#' distributions and get fitted values, confidence intervals, 
#' posterior predictive draws, and posterior predictive intervals
#' 
#' here, y, X are from the probability sample
#' and mu_beta, ..., b are hyperparameters
nps_mc <- function(size, y, X, mu_beta, k0, V, a, b, typeIerror = 0.05) {
    n <- nrow(X)
    p <- ncol(X)

    XTX <- crossprod(X)
    inv_scale_mat <- solve(XTX + solve(k0 * V))
    
    beta_hat <- solve(XTX) %*% crossprod(X, y) # eqn 7
    W <- inv_scale_mat %*% XTX # eqn 6
    post_mn <- W %*% beta_hat + (diag(1, p) - W) %*% mu_beta # eqn 5

    RSS <- crossprod(y - X %*% beta_hat) # eqn 10
    bet_minus_mu <- beta_hat - mu_beta
    SS <- RSS + crossprod(bet_minus_mu, solve(solve(XTX) + k0 * V) %*% bet_minus_mu) # eqn 9
    SS <- SS[1] # convert to scalar data type
    post_var <- inv_scale_mat * (SS + 2 * b) / (n + 2 * a - 2) # eqn 8
    
    scale_mat <- post_var * (n + 2 * a - 2) / (n + 2 * a)
    d_f <- n + 2 * a

    post_mn <- as.vector(post_mn)

    post_bet <- mvtnorm::rmvt(size, scale_mat, d_f, post_mn, "shifted")
    bet_mean <- colMeans(post_bet)
    post_tau <- rgamma(size, n / 2 + a, SS /2 + b)
    tau_mean <- mean(post_tau)

    # post_bet is size x p
    # fitted_y is size x n
    # X is n * p
    fitted_y <- tcrossprod(post_bet, X) 
    fitted_y_point <- colMeans(fitted_y)
    fitted_y_intervals <- apply(fitted_y, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))

    post_pred_draws <- rnorm(size * n, as.vector(fitted_y), 
                                rep(sqrt(1 / post_tau), n))
    post_pred_draws <- matrix(post_pred_draws, nrow = size)
    post_pred_intervals <- apply(post_pred_draws, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))

    nps_list <- list(beta_draws = post_bet, 
                     beta_mean = bet_mean,
                     tau_draws = post_tau, 
                     tau_mean = tau_mean,
                     fitted_y_draws = fitted_y, 
                     fitted_y_means = fitted_y_point,
                     fitted_y_lower = fitted_y_intervals[1,], 
                     fitted_y_upper = fitted_y_intervals[2,],
                     post_pred_draws = post_pred_draws, 
                     post_pred_lower = post_pred_intervals[1,], 
                     post_pred_upper = post_pred_intervals[2,], 
                     typeIerror = typeIerror)
    
    class(nps_list) <- "npsPost"
    return(nps_list)
}

# Make predictions and output mean response credible intervals as well 
# as posterior predictive draws and intervals
predict.npsPost <- function(object, newdata, typeIerror = NULL) {
    beta_draws <- object$beta_draws
    tau_draws <- object$tau_draws
    size <- nrow(beta_draws)
    n <- nrow(newdata)

    y_pred_draws <- tcrossprod(beta_draws, newdata)
    y_pred_mean <- colMeans(y_pred_draws)

    if(is.null(typeIerror)) typeIerror <- object$typeIerror
    mean_pred_intervals <- apply(y_pred_draws, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))

    post_pred_draws <- rnorm(size * n, as.vector(y_pred_draws), 
                                rep(sqrt(1 / tau_draws), n))
    post_pred_draws <- matrix(post_pred_draws, nrow = size)
    post_pred_intervals <- apply(post_pred_draws, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))

    nps_list <- list(beta_draws = beta_draws, 
                     beta_mean = colMeans(beta_draws),
                     tau_draws = tau_draws, 
                     tau_mean = mean(tau_draws),
                     fitted_y_draws = y_pred_draws, 
                     fitted_y_means = y_pred_mean,
                     fitted_y_lower = mean_pred_intervals[1,], 
                     fitted_y_upper = mean_pred_intervals[2,],
                     post_pred_draws = post_pred_draws, 
                     post_pred_lower = post_pred_intervals[1,], 
                     post_pred_upper = post_pred_intervals[2,], 
                     typeIerror = typeIerror)
    
    class(nps_list) <- "npsPost"
    return(nps_list)
}
#' Calculate Adaptive Power Prior Exponent for Logistic Regression
#'
#' Uses Hotelling's T^2 test to determine the power prior exponent `a`
#' based on the agreement between PS and NPS coefficient estimates.
#' Following Salvatore et al. (2024) approach.
#'
#' @param y_ps Binary outcome vector from probability sample
#' @param X_ps Design matrix from probability sample (without intercept)
#' @param y_nps Binary outcome vector from nonprobability sample
#' @param X_nps Design matrix from nonprobability sample (without intercept)
#' @param verbose If TRUE, prints diagnostic information
#'
#' @return A list with:
#'   - a: The adaptive power prior exponent (p-value from Hotelling test)
#'   - t2: Hotelling's T^2 statistic
#'   - f_stat: F-statistic
#'   - p_value: p-value (same as `a`)
#'   - beta_ps: PS coefficient estimates
#'   - beta_nps: NPS coefficient estimates
#'   - diff: Difference between PS and NPS coefficients
#'
#' @details
#' Lower `a` values indicate greater disagreement between PS and NPS,
#' suggesting NPS should receive less weight. Higher values indicate
#' good agreement and allow more NPS influence.
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
