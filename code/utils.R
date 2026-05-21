# Draw from MVN by specifying a covariance matrix
rmvnorm = function(n, mean, covar)
{
  k <- length(mean)
  stopifnot(k == nrow(covar) && k == ncol(covar))
  Z <- matrix(rnorm(n*k), k, n)
  A <- t(chol(covar))
  out <- A %*% Z + mean

  if(n==1){ return(as.vector(out)) } else { return(out) }
}

# Draw from MVN by specifying a precision matrix
rmvnorm_prec = function(n, mean, prec)
{
  k <- length(mean)
  stopifnot(k == nrow(prec) && k == ncol(prec))
  Z <- matrix(rnorm(n*k), k, n)

  # Note that Ainv %*% t(Ainv) is the Cholesky decomposition of the covariance
  # matrix solve(Omega)
  A <- chol(prec)
  Ainv <- backsolve(A, diag(1,k,k))
  out <- Ainv %*% Z + mean

  if(n==1){ return(as.vector(out)) } else { return(out) }
}


int_score <- function(alpha, truth, L, U){
  return(
    (U - L) + 2/alpha*(truth < L)*(L - truth) + 2/alpha*(truth > U)*(truth - U)
  )
}

# Estimate pseudo‐inclusion weights for a non‐probability sample
# via a logistic model distinguishing prob vs. non‐prob cases

estimate_ipw <- function(ps, nps, cov_formula, method) {
  stopifnot(method %in% c("ignorable", "beta_reg", "weighted"))

  ps_weights <- ps$weights #Need to double check whether scaling matters here
      
  if(method == "beta_reg") {
    filter_id <- which(ps_weights == 1)
    ps$weights[filter_id] <- 1.001
    beta_reg_out <- betareg::betareg(formula=as.formula(paste("1/weights",
                                                     deparse(cov_formula))),
                            data=ps)
    nps_prob_pred <- predict(beta_reg_out, newdata=nps)    
  } else if(method == "weighted") {
    combined <- bind_rows(
      ps  %>% mutate(Z = 0),
      nps %>% mutate(Z = 1)
    )

    combined <- combined %>%
        mutate(weights = ifelse(Z == 1, 1, weights))         
        
    fit_prop <- glm(
      formula = as.formula(paste0("Z ", deparse(cov_formula))),
      weights = weights, 
      family  = binomial(),
      data    = combined
    )

    nps_prob_pred  <- predict(fit_prop, newdata = nps, type = "response")
  }

  # NOTE: "ignorable" method is accepted but not implemented. Only "beta_reg" and
  # "weighted" are currently used. If "ignorable" is called, nps_prob_pred will
  # be undefined and will error on the next line.
  nps_weights <- 1/nps_prob_pred

  # Normalize 
  n_ps <- nrow(ps)
  n_nps <- nrow(nps)
  C_s <- n_ps/(n_nps + n_ps)
  C_s_star <- n_nps/(n_nps + n_ps)

  nps_weights <- nps_weights * (C_s_star * sum(nps_weights)/sum(ps_weights))
  ps_weights <- ps_weights * C_s

  return(list(nps_ipw = nps_weights,
              ps_ipw = ps_weights))
}

# Horvitz-Thompson (HT) Direct Estimator
#
# Calculates direct estimates by PUMA using survey-weighted Horvitz-Thompson estimator.
# Used for IPW-based direct estimation where weights combine PS and NPS samples.
#
# @param df Data frame with columns: response variable, PUMA (area), weights
# @param model_name String to label the model in output
# @param response_var Name of response variable
# @param alpha Significance level for confidence intervals
# @return Data frame with PUMA-level estimates: PUMA, point_est, lower_CI, upper_CI, model
HT <- function(df, model_name, response_var = NULL, alpha = 0.05) {
  # If response_var not provided, try to get from parent environment
  if (is.null(response_var)) {
    response_var <- get("response_var", envir = parent.frame())
  }

  samp.design <- survey::svydesign(ids = ~1, weights = ~weights, data = df)

  # Calculate direct estimates by PUMA with standard errors
  direst <- survey::svyby(
    as.formula(paste0("~", response_var)),
    ~PUMA,
    samp.design,
    survey::svymean,
    na.rm = TRUE,
    vartype = "se",
    keep.names = FALSE
  ) %>%
    dplyr::arrange(PUMA) %>%
    dplyr::transmute(
      PUMA,
      point_est = .data[[response_var]],
      lower_CI = point_est + qnorm(alpha / 2) * se,
      upper_CI = point_est + qnorm(1 - alpha / 2) * se,
      model = model_name
    )

  return(direst)
}
