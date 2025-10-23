#' NPS Prior Method

#' Important: marginal posterior of beta is multivariate t with 
#' location `post_mn`, scale matrix `post_var * (n + 2 * a - 2) / (n + 2 * a)`, 
#' and degrees of freedom = `n + 2 * a`
#' 
#' Generate Monte Carlo samples from the marginal posterior 
#' distributions and get fitted values and confidence intervals
#' for area-level predictions
 
#' Helper function to construct candidate prior distributions
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
        crossprod() %>% 
        as.numeric() / (n_NP - p)
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
        crossprod() %>% 
        as.numeric()
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
        name = "Noninformative",
        mu_beta = rep(0, p),
        k0 = n,
        V = I_p
    )

    conj <- list(
        name = "Conjugate",
        mu_beta = beta_hat_NP,
        k0 = if(pval < 0.05) 1 / log(n_NP) else 1 / n_NP,
        V = I_p
    )

    Zellner <- list(
        name = "Zellner",
        mu_beta = beta_hat_NP,
        k0 = if(pval < 0.05) n_NP^2 else 1,
        V = X_NP_cross_inv
    )

    Vd <- pmax(abs(diff_vec), sqrt(beta_hat_NP_var)) %>% 
        as.numeric() %>% 
        diag()
    conj_dist <- list(
        name = "Conjugate-distance",
        mu_beta = beta_hat_NP,
        k0 = 1 / log(n_NP),
        V = Vd^2
    )

    Zellner_dist <- list(
        name = "Zellner-distance",
        mu_beta = beta_hat_NP,
        k0 = n_NP,
        V = Vd %*% X_NP_cross_inv %*% Vd
    )

    # New options
    unit_info_NP <- list(
        name = "Unit Info NP",
        mu_beta = beta_hat_NP,
        k0 = n_NP,
        V = X_NP_cross_inv
    )

    unit_info_P <- list(
        name = "Unit Info P",
        mu_beta = beta_hat_NP,
        k0 = n,
        V = X_NP_cross_inv
    )

    ratio_info <- list(
        name = "Ratio Info",
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
        X = X,
        beta_hat_P = beta_hat_P, 
        X_cross_P = X_cross,
        X_cross_inv_P = X_cross_inv, 
        RSS_P = beta_hat_P_RSS
    )

    return(res)
}

#' Workhorse/internal function which does the Monte Carlo 
#' approximation for each candidate prior
nps_mc <- function(
    niter, 
    prior, 
    n,
    X,
    beta_hat_P, 
    X_cross_P,
    X_cross_inv_P, 
    RSS_P,
    PUMA_levels,
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
    SS <- RSS_P + crossprod(bet_minus_mu, chol2inv(chol(X_cross_inv_P + k0 * V)) %*% bet_minus_mu) # eqn 9
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
    # X is n * p
    # mean_respose is niter x n
    mean_response <- tcrossprod(post_bet, X) 
    post_pred_draws <- rnorm(niter * n, as.vector(mean_response), 
                                rep(sqrt(1 / post_tau), n))
    post_pred_draws <- matrix(post_pred_draws, nrow = niter)

    # aggregate by area
    by_PUMA_means <- sapply(PUMA_levels, function(cur_puma) {
      rowMeans(post_pred_draws[, which(PUMA == cur_puma), drop = FALSE])
    })

    point_est <- colMeans(by_PUMA_means)
    CI <- colQuantiles(by_PUMA_means, probs = c(0.025, 0.975))

    res_df <- data.frame(
        PUMA = PUMA_levels, 
        point_est = point_est, 
        lower_CI = CI[,1], 
        upper_CI = CI[,2], 
        submodel = prior$name
    )

    return(res_df)
}

#' This is the function the user should call
#' 
#' It returns a data frame with the area-level 
#' point and interval estimates 
#' for each candidate prior
nps_fit <- function(
    niter, 
    y_lik, 
    X_lik, 
    y_prior, 
    X_prior, 
    PUMA_levels,
    typeIerror = 0.05, 
    a = 0, 
    b = 0
) {
    prior_info <- make_priors(y_lik, X_lik, y_prior, X_prior)

    prior_list <- prior_info$prior_list

    res_df <- lapply(prior_list, function(prior) {
        nps_mc(
            niter, 
            prior, 
            prior_info$n,
            prior_info$X,
            prior_info$beta_hat_P, 
            prior_info$X_cross_P,
            prior_info$X_cross_inv_P, 
            prior_info$RSS_P, 
            PUMA_levels,
            typeIerror,
            a,
            b
        )
    }) %>% do.call(rbind, .) %>% 
            mutate(model = "Wis", .before = "submodel")

    return(res_df)
}