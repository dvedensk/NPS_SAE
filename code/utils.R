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

  ps_weights <- ps$PWGTP #Need to double check whether scaling matters here
      
  if(method == "beta_reg") {
    filter_id <- which(ps_weights == 1)
    ps$PWGTP[filter_id] <- 1.001
    beta_reg_out <- betareg::betareg(formula=as.formula(paste("1/PWGTP",
                                                     deparse(X_formula))),
                            data=ps)
    nps_prob_pred <- predict(beta_reg_out, newdata=nps)    
  } else if(method == "weighted") {
    combined <- bind_rows(
      ps  %>% mutate(Z = 0),
      nps %>% mutate(Z = 1)
    )

    combined <- combined %>%
        mutate(PWGTP = ifelse(Z == 1, 1, PWGTP))         
        
    fit_prop <- glm(
      formula = as.formula(paste0("Z ", deparse(cov_formula))),
      weights = PWGTP, 
      family  = binomial(),
      data    = combined
    )

    nps_prob_pred  <- predict(fit_prop, newdata = nps, type = "response")
  }

  # FIXME: if method == ignorable (used in main), nps_prob_pred is never defined
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
