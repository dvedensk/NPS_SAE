mixed_dist_mcmc <- function(
    X, 
    wts,
    PUMA_levels,
    prior_prec,
    beta_NP_hat, 
    kappas,
    prior_name
) {
    n <- nrow(X)
    p <- ncol(X)
    bet <- rep(0, p)
    y_pred <- matrix(NA, niter, n)

    for(s in 1:niter) {
        # Polya-Gamma draws
        omega <- rpg(n, wts, as.numeric(X %*% bet))

        # Calculate posterior covariance and mean of coefficients
        post_prec <- crossprod(X, diag(omega)) %*% X + prior_prec
        post_cov <- post_prec %>% 
            chol() %>% 
            chol2inv()

        post_mn <- as.numeric(post_cov %*% (crossprod(X, kappas) + prior_prec %*% beta_NP_hat))

        # Draw coefficients
        bet <- rmvnorm(1, post_mn, post_cov)[1,]

        # Draw from the ppd
        y_pred[s,] <- rbinom(n, 1, 1 / (1 + exp(-as.numeric(X %*% bet))))
    }

    # Aggregate by area
    by_PUMA_means <- sapply(PUMA_levels, function(cur_puma) {
      rowMeans(post_pred_draws[, which(PUMA == cur_puma), drop = FALSE])
    })

    # TODO: test equalite
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
        submodel = prior_name
    )

    return(res_df)
}

power_prior_mcmc <- function(
    X,
    X_NP,
    wts, 
    a,
    PUMA_levels,
    kappas_all, 
    niter
) {
    n <- nrow(X)
    n_NP <- nrow(X_NP)
    p <- ncol(X)

    X_all <- rbind(X, X_all)

    bet <- rep(0, p)
    y_pred <- matrix(NA, niter, n)

    for(s in 2:niters) {
        # Latent variance draws
        tau2 <- 1 / rgamma(p, 2, (18.75 + bet^2) / 2)

        # PS Polya Gamma Draws
        omega <- rpg(n, wts, as.numeric(X %*% bet))
        # NPS Polya Gamma Draws
        omega_NP <- rpg(n_NP, a, as.numeric(X_NP %*% bet))

        # Calculate posterior covariance and mean of regression coefficients
        post_prec <- crossprod(X_all, diag(c(omega, omega_NP))) %*% X_all + 
            diag(1 / tau2)

        post_cov <- post_prec %>%
            chol() %>% 
            chol2inv()

        post_mn <- as.numeric(post_cov %*% crossprod(X_all, kappas_all))

        # Draw standard normals, scale and shift to the posterior distribution
        bet <- rmvnorm(1, post_mn, post_cov)[1,]

        # Draw from the ppd
        y_pred[s,] <- rbinom(n, 1, 1 / (1 + exp(-as.numeric(X %*% bet))))
    }

    # Aggregate by area
    by_PUMA_means <- sapply(PUMA_levels, function(cur_puma) {
      rowMeans(post_pred_draws[, which(PUMA == cur_puma), drop = FALSE])
    })

    # TODO: test equalite
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
        submodel = "Power Prior"
    )

    return(res_df)
}

# MCMC outer function
nps_prior_mcmc <- function(
    y, 
    X, 
    y_NP, 
    X_NP, 
    niter, 
    PUMA_levels,
    wts = NULL, 
    typeIerr = 0.05
) {
    browser()
    n <- length(y)
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

    # NPS
    glm_fit_NP <- glm(y_NP ~ X_NP - 1, binomial())
    beta_NP_hat <- coef(glm_fit_NP) %>% 
        as.numeric()
    
    diff_vec <- beta_hat - beta_NP_hat
    sq_dist <- diff_vec^2
    kappas <- wts * (y - 0.5)

    # Mixed-distance prior
    prior_prec <- diag(1 / sq_dist)

    md_res <- mixed_dist_mcmc(
        X, 
        wts,
        PUMA_levels,
        prior_prec,
        beta_NP_hat, 
        kappas, 
        "Mixed-Distance"
    )

    # for mixed-distance-log priors, we estimate the variance of the MLE of each coefficient
    # based on the NPS 
    NP_vars <- glm_fit_NP %>% 
        vcov() %>% 
        diag() %>% 
        as.numeric()

    maxes <- pmax(sq_dist, NP_vars) 

    # Mixed-distance-log
    prior_vars <- maxes / log(n_NP)
    prior_prec <- diag(1 / prior_vars)

    mdl_res <- mixed_dist_mcmc(
        X, 
        wts,
        PUMA_levels,
        prior_prec,
        beta_NP_hat, 
        kappas, 
        "Mixed-Distance-log"
    )

    # Mixed_distance-log10
    prior_vars <- maxes / log10(n_NP)
    prior_prec <- diag(1 / prior_vars)

    mdl_ten_res <- mixed_dist_mcmc(
        X, 
        wts,
        PUMA_levels,
        prior_prec,
        beta_NP_hat, 
        kappas, 
        "Mixed-Distance-log10"
    )

    # Power Prior
    # Set power prior exponent to p-value corresponding to Hotelling T^2 test
    cov_beta_hat <- vcov(glm_fit)
    t2 <- crossprod(diff_vec, (cov_beta_hat %>% chol() %>% chol2inv()) %*% diff_vec)
    f = p * (n - 1) / (n - p)
    Fstat = t2 / f
    a <- pf(Fstat, p, n - p, lower.tail = F)

    kappas_all <- c(kappas, a * (y_NP - 1/2))

    pp_res <- power_prior_mcmc(
        X,
        X_NP,
        wts, 
        a,
        PUMA_levels,
        kappas_all, 
        niter
    )

    res_df <- rbind(
        md_res, 
        mdl_res, 
        mdl_ten_res,
        pp_res
    )

    res_df <- res_df %>% 
        mutate(model = "NPS Prior", .before = "submodel")

    return(res_df)
}