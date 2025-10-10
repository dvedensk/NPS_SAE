library(BayesLogit)
library(Matrix)
library(LaplacesDemon)

# This function fits a survey-weighted, binary unit-level model via Gibbs sampling
# X is the input covariate matrix, one row per sample unit
# Psi is the matrix of spatial basis functions, one row per sample unit
# y is a vector of binary survey responses
# sigma2_beta is the prior variance for Beta
# iter is the number of iterations
# burn is the length of burn in
# weights is the vector of survey weights
fit_bulm <- function(X, Psi, y, sigma2_beta=1000, iter=1000, burn=500, weights=NULL){
  p <- ncol(X)
  r <- ncol(Psi)
  n <- length(y)
  if(is.null(weights)) weights <- rep(1, n)
  b_inv <- (1/sigma2_beta) * Diagonal(p)
  kappa <- weights * (y - 0.5)
  w <- rep(1, n)
  beta <- rep(0, p)
  eta <- rep(0, r)
  sigma2_eta <- 1
  beta_out <- matrix(NA, nrow=iter, ncol=p)
  eta_out <- matrix(NA, nrow=iter, ncol=r)
  sigma2_eta_out <- rep(NA, iter)

  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## Sample fixed effects
    prec_beta <- solve(t(X) %*% Diagonal(length(w), w) %*% X + b_inv)
    mean_beta <- t(X) %*% Diagonal(length(w), w) %*% (kappa/w - Psi %*% eta )
    beta <- beta_out[i, ]  <- as.numeric(rmvnorm(1, 
						mean=prec_beta %*% mean_beta, 
						sigma=as.matrix(prec_beta)))
    
    ## Sample random effects with sum-to-zero constraint
    Einv <- (1/sigma2_eta) * Diagonal(r)
    prec_eta <- solve(t(Psi) %*% Diagonal(length(w), w) %*% Psi + Einv)
    mean_eta <- t(Psi) %*% Diagonal(length(w),w) %*% (kappa/w - X %*% beta )
    eT <- as.numeric(rmvnorm(1, 
			     mean=prec_eta %*% mean_eta, 
			     sigma=as.matrix(prec_eta)))
    eta <- eta_out[i,]  <- eT - mean(eT)
    
    ## Sample RE variance
    sigma2_eta <- sigma2_eta_out[i] <- 1/rgamma(1,
                                     shape=0.5 + r/2,
                                     (0.5 + 0.5 * t(eta) %*% eta))
    
    ## Sample latent PG variables
    w <- rpg(n, weights, as.numeric(X%*%beta + Psi%*%eta))
    setTxtProgressBar(pb, i)
  }
  return(list(beta=beta_out[-c(1:burn),],
              eta=eta_out[-c(1:burn),],
              sigma2_eta=sigma2_eta_out[-c(1:burn)]))
}

post_preds <- function(grouped_pop_df, beta, eta, alpha, X_formula, Psi_formula) {
  X_pop <- model.matrix(X_formula , data=grouped_pop_df)
  Psi_pop <- model.matrix(Psi_formula, data=grouped_pop_df)
  pop_size <- nrow(grouped_pop_df)
  probs <- plogis(X_pop %*% t(beta) + Psi_pop %*% t(eta))
  preds <- rbinom(n=pop_size*nrow(beta), size=grouped_pop_df$n, prob=c(probs))
  preds <- matrix(preds, nrow=pop_size)
  agg_df <- data.frame(grouped_pop_df[,c("PUMA", "n")], MCMC=preds) %>%
      group_by(PUMA) %>%
      summarize_all(sum)
  #rescale to proportions:
  agg_df[,-c(1,2)] <- agg_df[,-c(1,2)]/agg_df$n 
  summary_df <- agg_df %>% rowwise(PUMA) %>%
      summarize(point_est=mean(c_across(starts_with("MCMC"))),
                lower_CI=quantile(c_across(starts_with("MCMC")), alpha/2),
                upper_CI=quantile(c_across(starts_with("MCMC")), 1 - alpha/2))
  return(summary_df)
}

bulm_results <- function(grouped_pop_df, alpha, X, Psi, y, sigma2_beta=1000,
                         X_formula, Psi_formula, iter=1000, burn=500,
                         weights=NULL, summaries=TRUE) {
  coeffs <- fit_bulm(X, Psi, y, sigma2_beta, iter, burn, weights)
  if (summaries){
    return(post_preds(grouped_pop_df, coeffs$beta, coeffs$eta,
                      alpha, X_formula, Psi_formula))
  } else{
    return(coeffs)
  }
}
