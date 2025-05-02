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