#' NPS Prior Method

#' Do a Bayesian linear regression with several variations on a 
#' Zellner's g prior on the regression coefficients given the precision. 
#' This form uses the probability sample in the likelihood and the 
#' nonprobability sample in the g prior. There is also a Gamma prior 
#' on the precision parameter.
#' 
#' We draw samples from the marginal posteriors 
#' of the regression coefficients and the precision. These distributions 
#' are specified in the paper. For reference here, the coefficients have 
#' a multivariate t posterior with  location `post_mn`, scale matrix 
#' `scale_mat`, and degrees of freedom = `n + 2a`. 
#' The precision has a Gamma posterior with shape `a + n / 2` and rate 
#' `b + SS / 2`.
#' 
#' Generate Monte Carlo samples from the marginal posterior 
#' distributions and get fitted values and confidence intervals
#' for area-level predictions
 
#' Helper function to construct hyperparameters for the prior distributions.
#' A lot of the prior distributions make the same calculations, so constructing 
#' all priors in one function prevents redundant computation.
make_priors <- function(y, X, y_NP, X_NP, wts) {
    n <- nrow(X)
    p <- ncol(X)
    n_NP <- nrow(X_NP)

    # NP sample: regression coefficients, MSE, cov(beta_hat), var(beta_hat_j) for all j = 1, ..., p
    X_NP_cross <- crossprod(X_NP)
    X_NP_cross_inv <- X_NP_cross %>% 
        chol() %>% 
        chol2inv()

    beta_hat_NP <- X_NP_cross_inv %*% crossprod(X_NP, y_NP)
    beta_hat_NP_MSE <- (y_NP - X_NP %*% beta_hat_NP) %>% 
        crossprod() %>% 
        as.numeric() / (n_NP - p)
    vcovm_beta_hat_NP <- beta_hat_NP_MSE * X_NP_cross_inv
    beta_hat_NP_var <- diag(vcovm_beta_hat_NP)

    # P sample: egression coefficients, RSS (residual SS), MSE, cov(beta_hat)

    # Pseudolikelihood
    Wt_mat_sqrt <- diag(sqrt(wts))
    y_wt <- Wt_mat_sqrt %*% y
    X_wt <- Wt_mat_sqrt %*% X 
    
    X_cross <- crossprod(X_wt)
    X_cross_inv <- X_cross %>% 
        chol() %>% 
        chol2inv()

    beta_hat_P <- X_cross_inv %*% crossprod(X_wt, y_wt)
    beta_hat_P_RSS <- (y - X %*% beta_hat_P) %>% 
        crossprod() %>% 
        as.numeric()
    beta_hat_P_MSE <- beta_hat_P_RSS / (n - p)

    # Calculate the covariance matrix of the regression coefficients. The extra complexity is due to using a pseudolikelihood; 
    # it's like using a WLS estimate for beta_hat_P but assuming Var(y) = sigma^2 I_n so the extra terms don't cancel.
    # If the weights are all = 1, then this reduces to the typical OLS estimate of MSE * (X'X)^{-1}
    vcovm_beta_hat_P <- beta_hat_P_MSE * X_cross_inv %*% crossprod(X, Wt_mat_sqrt^4 %*% X) %*% X_cross_inv

    # Hotelling T^2 test
    diff_vec <- beta_hat_P - beta_hat_NP

    t2 <- crossprod(diff_vec, (vcovm_beta_hat_P %>% chol() %>% chol2inv()) %*% diff_vec)
    f = p * (n - 1) / (n - p)
    Fstat = t2 / f

    pval = pf(Fstat, p, n - p, lower.tail = F)

    # Construct priors
    I_p <- diag(1, p)

    # Options listed in the paper
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

    # If we make an R package, it is useful to specify a class template.
    # Else, we can remove this chunk.
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

    # Return the list of hyperparameters as well as other 
    # statistics which will be used when calculating posterior parameters.
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

#' Workhorse/internal function which does the Monte Carlo 
#' draws for a single candidate prior
nps_mc <- function(
    niter, 
    prior, 
    X,
    n,
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

    # Load prior parameters
    mu_beta <- prior$mu_beta 
    k0 <- prior$k0
    V <- prior$V

    # Calculate posterior parameters
    inv_scale_mat <- (X_cross_P + chol2inv(chol(k0 * V))) %>% 
        chol() %>% 
        chol2inv()

    W <- inv_scale_mat %*% X_cross_P # eqn 6
    post_mn <- W %*% beta_hat_P + (diag(1, p) - W) %*% mu_beta # eqn 5

    bet_minus_mu <- beta_hat_P - mu_beta
    SS <- RSS_P + crossprod(bet_minus_mu, chol2inv(chol(X_cross_inv_P + k0 * V)) %*% bet_minus_mu) # eqn 9
    SS <- as.numeric(SS) # convert to scalar data type
    d_f <- n + 2 * a 
    post_var <- inv_scale_mat * (SS + 2 * b) / (d_f - 2) # eqn 8
    
    scale_mat <- post_var * (d_f - 2) / (d_f)

    post_mn <- as.vector(post_mn)

    # Draw from posteriors
    post_bet <- rmvt(niter, scale_mat, d_f, post_mn, "shifted")
    post_tau <- rgamma(niter, a + n / 2, b + SS/2)

    # Draw unit-level predictions from posterior predictive distribution
    mean_response <- tcrossprod(post_bet, X) 
    post_pred_draws <- rnorm(
        niter * n, 
        as.vector(mean_response), 
        rep(sqrt(1 / post_tau), n)
    ) %>% 
        matrix(nrow = niter)

    # Aggregate by area
    by_PUMA_means <- sapply(PUMA_levels, function(cur_puma) {
      rowMeans(post_pred_draws[, which(PUMA == cur_puma), drop = FALSE])
    })

    # This is an untested chunk which might be a lot faster than using sapply
    # by_PUMA_means2 <- post_pred_draws %>% 
    #     t() %>% 
    #     rowsum(group = PUMA, reorder = FALSE) %>% 
    #     t()
    # 
    # by_PUMA_means2 <- py_PUMA_means2 / as.numeric(table(PUMA))

    # Calculate area-level point and interval estimates
    point_est <- colMeans(by_PUMA_means)
    CI <- colQuantiles(by_PUMA_means, probs = c(typeIerror / 2, 1 - typeIerror / 2))

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
    wts = NULL,
    typeIerror = 0.05, 
    a = 0, 
    b = 0
) {
    n <- length(y_lik)
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

    prior_info <- make_priors(y_lik, X_lik, y_prior, X_prior, wts)

    prior_list <- prior_info$prior_list

    res_df <- lapply(prior_list, function(prior) {
        nps_mc(
            niter, 
            prior, 
            X_lik,
            n,
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
