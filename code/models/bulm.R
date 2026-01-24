library(BayesLogit)
library(Matrix)

# This function fits a survey-weighted, binary unit-level model via Gibbs sampling
# X is the input covariate matrix, one row per sample unit
# Psi is the matrix of spatial basis functions, one row per sample unit
# y is a vector of binary survey responses
# sigma2_beta is the prior variance for Beta
# iter is the number of iterations
# burn is the length of burn in
# weights is the vector of survey weights
fit_bulm <- function(X, Psi, y, sigma2_beta=1000, iter=1000, burn=500,
                     weights=NULL, display_progress=FALSE){
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

  if(display_progress) { pb <- txtProgressBar(min=0, max=iter, style=3) }
  for(i in 1:iter){

    ## Sample fixed effects
    prec_beta <- crossprod(X * sqrt(w)) + b_inv
    mean_beta <- t(X) %*% Diagonal(length(w), w) %*% (kappa/w - Psi %*% eta )
    beta <- beta_out[i, ]  <- rmvnorm_prec(n=1, 
					   mean=qr.solve(prec_beta, mean_beta),
					   prec=prec_beta)
    
    ## Sample random effects with sum-to-zero constraint
    prec_eta <- crossprod(Psi * sqrt(w)) + Diagonal(r)/sigma2_eta
    mean_eta <- t(Psi) %*% Diagonal(length(w),w) %*% (kappa/w - X %*% beta )
    eT <- rmvnorm_prec(n=1, 
	   	       mean=qr.solve(prec_eta, mean_eta),
		       prec=prec_eta)
    eta <- eta_out[i,]  <- eT - mean(eT)
    
    ## Sample RE variance
    sigma2_eta <- sigma2_eta_out[i] <- 1/rgamma(1,
                                                shape=0.5 + r/2,
                                                (0.5 + 0.5 * t(eta) %*% eta))
    
    ## Sample latent PG variables
    w <- rpg(n, weights, as.numeric(X%*%beta + Psi%*%eta))
    if(display_progress) { setTxtProgressBar(pb, i) }
  }
  return(list(beta=beta_out[-c(1:burn),],
              eta=eta_out[-c(1:burn),],
              sigma2_eta=sigma2_eta_out[-c(1:burn)]))
}

post_preds <- function(grouped_pop_df, beta, eta, alpha, X_formula, Psi_formula, stan=FALSE) {
  if(stan) {
    beta <- posterior::as_draws_matrix(beta)
    eta <- posterior::as_draws_matrix(eta)
  }

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

generate_bulm_results <- function(grouped_pop_df, alpha, X, Psi, y, sigma2_beta=1000,
                                  X_formula, Psi_formula, iter=1000, burn=500, weights=NULL,
                                  display_progress=FALSE) {
    
  coeffs <- fit_bulm(X, Psi, y, sigma2_beta, iter, burn,
                     weights, display_progress=display_progress)
  return(list(summaries = post_preds(grouped_pop_df, coeffs$beta, coeffs$eta,
                                     alpha, X_formula, Psi_formula),
              chains=coeffs))
}

get_stan_summaries <- function(y, X, Psi, weights, sigma2_beta=3, n_chains,
                               mcmc_burn, mcmc_iter, grouped_pop_df, alpha,
                               X_formula, Psi_formula, threads_per_chain=4,
                               seed = NULL) {
    
  bulm_stan_dat <- list(r=ncol(Psi), nn=length(y), p=ncol(X),
                        Y=y, weights=weights,
                        puma=apply(Psi, 1, which.max),
                        X=X, sigma2_beta=3)

  sample_args <- list(
    data = bulm_stan_dat,
    chains = n_chains,
    parallel_chains = n_chains,
    iter_warmup = mcmc_burn,
    iter_sampling = mcmc_iter,
    threads_per_chain = threads_per_chain
  )
  if (!is.null(seed)) sample_args$seed <- seed

  bulm_stan_out <- do.call(stan_bulm_mod$sample, sample_args)
   
  bulm_stan_pp <- post_preds(grouped_pop_df=acs_pop_grouped,
                             beta=bulm_stan_out$draws("beta"),
                             eta=bulm_stan_out$draws("eta"),
                             alpha=alpha,
                             X_formula=X_formula,
                             Psi_formula=Psi_formula,
                             stan=TRUE)

  return(bulm_stan_pp)
}
