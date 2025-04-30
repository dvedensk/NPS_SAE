# NPS Prior Method
# here, y, X are from the probability sample
# and mu_beta, ..., b are hyperparameters
nps_post <- function(y, X, mu_beta, k0, V, a, b) {
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

    post_mn <- as.vector(post_mn)

    fitted_y <- as.vector(X %*% post_mn)
    return(list(mn = post_mn, vr = post_var, fitted_y = fitted_y))
}