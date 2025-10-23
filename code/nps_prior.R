make_priors <- function(y, X, y_NP, X_NP) {
    n <- nrow(X)
    p <- ncol(X)
    n_NP <- nrow(X_NP)

    # Linear model fit for NP sample
    X_NP_cross <- crossprod(X_NP)
    X_NP_cross_inv <- X_NP_cross %>% 
        chol() %>% 
        chol2inv()
    
    # regression coefficients, MSE, cov(beta_hat), var(beta_hat_j) for all j = 1, ..., p
    beta_hat_NP <- X_NP_cross_inv %*% crossprod(X_NP, y_NP)
    beta_hat_NP_MSE <- (y_NP - X_NP %*% beta_hat_NP) %>% 
        crossprod()[1] / (n_NP - p)
    vcovm_beta_hat_NP <- beta_hat_NP_MSE * X_NP_cross_inv
    beta_hat_NP_var <- diag(vcovm_beta_hat_NP)

    # Linear model fit for P sample
    X_cross <- crossprod(X)
    X_cross_inv <- X_cross %>% 
        chol() %>% 
        chol2inv()

    # regression coefficients, RSS (residual SS), MSE, cov(beta_hat)
    beta_hat_P <- X_cross_inv %*% crossprod(X, y)
    beta_hat_P_RSS <- (y - X %*% beta_hat_NP) %>% 
        crossprod()[1]
    beta_hat_P_MSE <- beta_hat_P_RSS / (n - p)
    vcovm_beta_hat_P <- beta_hat_P_MSE * X_cross_inv

    # Hotelling T^2 test
    diff_vec <- beta_hat_P - beta_hat_NP

    t2 <- crossprod(diff_vec, (vcovm_beta_hat_P %>% chol() %>% chol2inv()) %*% diff_vec)
    f = p * (n - 1) / (n - p)
    Fstat = t2 / f

    pval = pf(Fstat, p, n - p, lower.tail = F)

    # Construct priors
    I_p <- diag(1, p)

    noninf <- list(
        mu_beta = rep(0, p),
        k0 = n,
        V = I_p
    )

    conj <- list(
        mu_beta = beta_hat_NP,
        k0 = if(test_res$pval < 0.05) 1 / log(n_NP) else 1 / n_NP,
        V = I_p
    )

    Zellner <- list(
        mu_beta = beta_hat_NP,
        k0 = if(test_res$pval < 0.05) n_NP^2 else 1,
        V = X_NP_cross_inv
    )

    Vd <- diag(max(abs(diff_vec), sqrt(beta_hat_NP_var)))
    conj_dist <- list(
        mu_beta = beta_hat_NP,
        k0 = 1 / log(n_NP),
        V = Vd^2
    )

    Zellner_dist <- list(
        mu_beta = beta_hat_NP,
        k0 = n_NP,
        V = Vd %*% X_NP_cross_inv %*%% Vd
    )

    # New options
    unit_info_NP <- list(
        mu_beta = beta_hat_NP,
        k0 = n_NP,
        V = X_NP_cross_inv
    )

    unit_info_P <- list(
        mu_beta = beta_hat_NP,
        k0 = n,
        V = X_NP_cross_inv
    )

    ratio_info <- list(
        mu_beta = beta_hat_NP,
        k0 = n_NP / n,
        V = X_NP_cross_inv
    )

    class(noninf) <- 
        class(conj) <- 
        class(Zellner) <- 
        class(conj_dist) <- 
        class(Zellner_dist) <- 
        class(unit_info_NP) <- 
        class(unit_info_P) <- 
        class(ratio_info) <-
        "npsPrior"
    
    prior_list <- list(
        noninf = noninf, 
        conj = conj, 
        Zellner = Zellner, 
        conj_dist = conj_dist, 
        Zellner_dist = Zellner_dist,
        unit_info_NP = unit_info_NP, 
        unit_info_P = unit_info_P, 
        ratio_info = ratio_info
    )

    res <- list(
        prior_list = prior_list,
        n = n,
        beta_hat_P = beta_hat_P, 
        X_cross_P = X_cross,
        X_cross_inv_P = X_cross_inv, 
        RSS_P = beta_hat_P_RSS
    )

    return(res)
}

nps_mc <- function(
    iter, 
    prior, 
    n,
    beta_hat_P, 
    X_cross_P,
    X_cross_inv_P, 
    RSS_P,
    typeIerror,
    a,
    b
) {
    p <- length(beta_hat_P)

    mu_beta <- prior$mu_beta 
    k0 <- prior$k0
    V <- prior$V

    inv_scale_mat <- chol2inv(chol((X_cross_P + chol2inv(chol(k0 * V)))))
    
    W <- inv_scale_mat %*% X_cross_P # eqn 6
    post_mn <- W %*% beta_hat_P + (diag(1, p) - W) %*% mu_beta # eqn 5

    bet_minus_mu <- beta_hat_P - mu_beta
    SS <- RSS + crossprod(bet_minus_mu, chol2inv(chol(X_cross_inv_P + k0 * V)) %*% bet_minus_mu) # eqn 9
    SS <- SS[1] # convert to scalar data type
    post_var <- inv_scale_mat * (SS + 2 * b) / (n + 2 * a - 2) # eqn 8
    
    scale_mat <- post_var * (n + 2 * a - 2) / (n + 2 * a)
    d_f <- n + 2 * a

    post_mn <- as.vector(post_mn)

    post_bet <- mvtnorm::rmvt(niter, scale_mat, d_f, post_mn, "shifted")
    bet_mean <- colMeans(post_bet)
    post_tau <- rgamma(niter, n / 2 + a, SS /2 + b)
    tau_mean <- mean(post_tau)

    # post_bet is niter x p
    # fitted_y is niter x n
    # X is n * p
    fitted_y <- tcrossprod(post_bet, X) 
    fitted_y_point <- colMeans(fitted_y)
    fitted_y_intervals <- apply(fitted_y, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))

    post_pred_draws <- rnorm(niter * n, as.vector(fitted_y), 
                                rep(sqrt(1 / post_tau), n))
    post_pred_draws <- matrix(post_pred_draws, nrow = niter)
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

nps_fit <- function(niter, y_lik, X_lik, y_prior, X_prior, typeIerror = 0.05, a = 0, b = 0) {
    prior_info <- make_priors(y_lik, X_lik, y_prior, X_prior)

    prior_list <- nps_prior_info$prior_list

    res <- lapply(prior_list, function(prior) {
        nps_mc(
            niter, 
            prior, 
            prior_info$n,
            prior_info$beta_hat_P, 
            prior_info$X_cross_P,
            prior_info$X_cross_inv_P, 
            prior_info$RSS_P, 
            typeIerror,
            a,
            b
        )

        # TODO: aggregate up to domain level
    })

    return(res)
}

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
#' 
#' 

# Outdated version
# nps_mc <- function(size, y, X, mu_beta, k0, V, a, b, typeIerror = 0.05) {
#     n <- nrow(X)
#     p <- ncol(X)
# 
#     XTX <- crossprod(X)
#     inv_scale_mat <- solve(XTX + solve(k0 * V))
#     
#     beta_hat <- solve(XTX) %*% crossprod(X, y) # eqn 7
#     W <- inv_scale_mat %*% XTX # eqn 6
#     post_mn <- W %*% beta_hat + (diag(1, p) - W) %*% mu_beta # eqn 5
# 
#     RSS <- crossprod(y - X %*% beta_hat) # eqn 10
#     bet_minus_mu <- beta_hat - mu_beta
#     SS <- RSS + crossprod(bet_minus_mu, solve(solve(XTX) + k0 * V) %*% bet_minus_mu) # eqn 9
#     SS <- SS[1] # convert to scalar data type
#     post_var <- inv_scale_mat * (SS + 2 * b) / (n + 2 * a - 2) # eqn 8
#     
#     scale_mat <- post_var * (n + 2 * a - 2) / (n + 2 * a)
#     d_f <- n + 2 * a
# 
#     post_mn <- as.vector(post_mn)
# 
#     post_bet <- mvtnorm::rmvt(size, scale_mat, d_f, post_mn, "shifted")
#     bet_mean <- colMeans(post_bet)
#     post_tau <- rgamma(size, n / 2 + a, SS /2 + b)
#     tau_mean <- mean(post_tau)
# 
#     # post_bet is size x p
#     # fitted_y is size x n
#     # X is n * p
#     fitted_y <- tcrossprod(post_bet, X) 
#     fitted_y_point <- colMeans(fitted_y)
#     fitted_y_intervals <- apply(fitted_y, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))
# 
#     post_pred_draws <- rnorm(size * n, as.vector(fitted_y), 
#                                 rep(sqrt(1 / post_tau), n))
#     post_pred_draws <- matrix(post_pred_draws, nrow = size)
#     post_pred_intervals <- apply(post_pred_draws, 2, function(y) quantile(y, probs = c(typeIerror / 2, 1 - typeIerror / 2)))
# 
#     nps_list <- list(beta_draws = post_bet, 
#                      beta_mean = bet_mean,
#                      tau_draws = post_tau, 
#                      tau_mean = tau_mean,
#                      fitted_y_draws = fitted_y, 
#                      fitted_y_means = fitted_y_point,
#                      fitted_y_lower = fitted_y_intervals[1,], 
#                      fitted_y_upper = fitted_y_intervals[2,],
#                      post_pred_draws = post_pred_draws, 
#                      post_pred_lower = post_pred_intervals[1,], 
#                      post_pred_upper = post_pred_intervals[2,], 
#                      typeIerror = typeIerror)
#     
#     class(nps_list) <- "npsPost"
#     return(nps_list)
# }

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