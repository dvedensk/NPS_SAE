library(readr)
library(dplyr)
library(tidyverse)
library(sampling)
library(mvtnorm)
library(survey)
library(purrr)
library(BayesLogit)
library(Matrix)
library(LaplacesDemon)
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr is required; install it via install.packages('cmdstanr').")
}
library(cmdstanr)
# Validate CmdStan is discoverable; cmdstanr 0.9.0 lacks cmdstan_available()
tryCatch(
  cmdstanr::cmdstan_path(),
  error = function(e) stop("CmdStan is not installed or toolchain is unavailable; run cmdstanr::install_cmdstan().")
)

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R")) # utils.R must define estimate_ipw()
source(file.path("code", "models", "bulm.R"))
source(file.path("code", "models", "VSW.R"))
source(file.path("code", "models", "IPW.R"))

source(file.path("code", "models", "mrp_all.R"))
# load .stan if available (allows running non-Stan paths without failure)
stan_bulm <- file.path("code", "models", "bulm.stan")
stan_si2 <- file.path("code", "models", "si2.stan")
stan_mrp_int2 <- file.path("code", "models", "mrp_int2.stan")
stan_models_available <- file.exists(stan_si2) && file.exists(stan_mrp_int2)
if (!stan_models_available) {
  warning("Stan model files not found; skipping Stan compilation and Stan-based models.")
  mod <- modINT <- NULL
} else {
  mod <- cmdstan_model(
    stan_si2,
    cpp_options = list(stan_threads = TRUE)
  )

  modINT <- cmdstan_model(
    stan_mrp_int2,
    cpp_options = list(stan_threads = TRUE)
  )

  stan_bulm_mod <- cmdstan_model(
    stan_bulm,
    cpp_options = list(stan_threads = TRUE)
  )

  Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())
}


set.seed(99)

# ============================================================
# CONFIGURATION: Response variable and sampling weights
# ============================================================
# Response variable to analyze (must exist in ACS_NPS_pop.csv)
response_var <- "PUBCOV" # Default: PUBCOV (binary). Can also use HICOV, WAGP, etc.
response_type <- "binary" # PUBCOV is 0/1 in processed population

# PS sampling weight configuration (for PUBCOV: use WAGP=0.05, PWGTP=-0.2)
PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)

# NPS sampling weight configuration (for PUBCOV: use PWGTP=0.3, AGEP=0.7)
NPS_weight_config <- list(PWGTP = 0.3, AGEP = 0.7)
# ============================================================

# load population file
# Note: AGEP is continuous (for sampling weights), AGEP_binned is categorical (for modeling)
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv")) %>%
  mutate(
    AGEP_binned = factor(AGEP_binned),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# population values to compare against
true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    response_true = mean(.data[[response_var]], na.rm = TRUE),
    .groups = "drop"
  )

acs_pop_grouped <- acs_pop %>%
  group_by(PUMA, AGEP_binned, RAC1P, SEX) %>%
  tally()

X_formula <- as.formula("~ AGEP_binned + RAC1P + SEX")
Psi_formula <- as.formula("~ -1 + PUMA")

mcmc_iter <- 1000
mcmc_burn <- 1000
n_chains <- 2
alpha <- .05 
N_pop <- nrow(acs_pop)
pop_mean <- mean(acs_pop[[response_var]], na.rm = TRUE)
cat(
  "\nPopulation summary:",
  "\n- N =", N_pop,
  "\n- Response mean (", response_var, ") =", round(pop_mean, 4),
  "\n- PS weights:", paste(names(PS_weight_config), PS_weight_config, sep = "=", collapse = ", "),
  "\n- NPS weights:", paste(names(NPS_weight_config), NPS_weight_config, sep = "=", collapse = ", "),
  "\n\n"
)

# take samples
Nsim <- 1
prob_samples <- list()
nonprob_samples <- list()
results <- list()
summary_df_VSW <- list()

for (sim in 1:Nsim) {
  print(sim)

  # 1. Draw probability and nonprobability samples
  ps_sample <- get_strat_PS(pop_df = acs_pop, samp_frac = .005, weight_config = PS_weight_config)
  nps_sample <- get_NPS(pop_df = acs_pop, samp_frac = .05, weight_config = NPS_weight_config, internet_only = FALSE)

  # save the samples
  prob_samples[[sim]] <- ps_sample
  nonprob_samples[[sim]] <- nps_sample

  # check the PS weights
  # FIXME: should this be true? It isn't. If it need not be true, we could add more explanation in the comment 
  # or remove this line. Do we also need a scaled version of the weights which sum to the population size?
  cat("PS weight Check: sum(ps_weights) =", sum(ps_sample$weights), "vs N_pop =", N_pop, "\n")

  # Calculate and print sample diagnostics (DDC and ESS)
  diagnostics <- calculate_sample_diagnostics(
    pop_df = acs_pop,
    ps_sample = ps_sample,
    nps_sample = nps_sample,
    response_var = response_var,
    print_results = TRUE,
    sim_num = sim
  )

  # 2. Scale weights for pseudolikelihood models
  ps_scale_weights <- length(ps_sample$idx) * ps_sample$weights / sum(ps_sample$weights)
  nps_scale_weights <- length(nps_sample$idx) * nps_sample$weights / sum(nps_sample$weights)

  # 3. Extract PS and NPS data frames
  ps <- acs_pop[ps_sample$idx, ] %>%
    mutate(weights = ps_sample$weights)
  nps <- acs_pop[nps_sample$idx, ]

  # 4. Build design matrices for PS and NPS separately
  X_ps <- model.matrix(X_formula, data = ps)
  Psi_ps <- model.matrix(Psi_formula, data = ps)
  y_ps <- ps[[response_var]]

  X_nps <- model.matrix(X_formula, data = nps)
  Psi_nps <- model.matrix(Psi_formula, data = nps)
  y_nps <- nps[[response_var]]

  X <- rbind(X_ps, X_nps)
  Psi <- rbind(Psi_ps, Psi_nps)
  y <- c(y_ps, y_nps)
  PUMAs <- c(ps$PUMA, nps$PUMA)

  # 5. Direct estimate on probability sample
  samp.design <- svydesign(ids = ~1, weights = ~weights, data = ps)

  # svyby extraction: use built-in se to avoid manual column juggling
  direst <- svyby(
    as.formula(paste0("~", response_var)),
    ~PUMA,
    samp.design,
    svymean,
    na.rm = TRUE,
    vartype = "se",
    keep.names = FALSE
  ) %>%
    arrange(PUMA) %>%
    transmute(
      PUMA,
      point_est = .data[[response_var]], # pulling the mean column from svyby
      lower_CI = point_est + qnorm(alpha / 2) * se,
      upper_CI = point_est + qnorm(1 - alpha / 2) * se,
      model = "direst"
    )

  # 6. Fit unit-level model on probability sample
  bulm_stan_dat <- list(r=ncol(Psi_ps), nn=length(y_ps), p=ncol(X_ps),
                        Y=y_ps, weights=ps_scale_weights,
                        puma=apply(Psi_ps, 1, which.max),
                        X=X_ps, sigma2_beta=3)

  bulm_stan_out <- stan_bulm_mod$sample(data=bulm_stan_dat, chains=n_chains,
                                        parallel_chains=n_chains, iter_warmup=mcmc_burn, 
                                        iter_sampling=mcmc_iter, threads_per_chain=4)
   
  bulm_ps_out <- post_preds(grouped_pop_df=acs_pop_grouped,
                            beta=bulm_stan_out$draws("beta"),
                            eta=bulm_stan_out$draws("eta"),
                            alpha=alpha,
                            X_formula=X_formula,
                            Psi_formula=Psi_formula,
                            stan=TRUE) 
  bulm_ps_out$model <- "bulm_ps_only"
  

  # 7. IPW
  # 7.A Estimate (two types of) IP weights for NPS
  weights_beta <- estimate_ipw(ps, nps, X_formula, "beta_reg")                            
  weights_uncond <- estimate_ipw(ps, nps, X_formula, "weighted")                                                                                                                        
  weights_beta <- c(weights_beta$ps_weights, weights_beta$nps_weights)                          
  weights_uncond <- c(weights_uncond$ps_weights, weights_uncond$nps_weights)    
  scale_weights_beta <- weights_beta/sum(weights_beta) * length(y) 
  scale_weights_uncond <- weights_uncond/sum(weights_uncond) * length(y) 

  #7.B Calculate Horvitz-Thompson (HT) style direct estimates using IPW
  beta_df <- data.frame(PUBCOV=y, PUMA=PUMAs, weights=weights_beta) 
  uncond_df <- data.frame(PUBCOV=y, PUMA=PUMAs, weights=weights_uncond)
  beta_HT <- HT(beta_df, "beta_HT")
  uncond_HT <- HT(uncond_df, "uncond_HT")
  
  beta_HT$model <- "IPW HT (beta)"
  uncond_HT$model <- "IPW HT (uncond)"

  #7.C Fit a unit-level model with IPW weights
  bulm_beta_out <- get_stan_summaries(y=y, X=X, Psi=Psi, weights=scale_weights_beta,
                                      n_chains=2, mcmc_burn=mcmc_burn,
                                      mcmc_iter=mcmc_iter, alpha=alpha,
                                      grouped_pop_df=acs_pop_grouped, 
                                      X_formula=X_formula, Psi_formula=Psi_formula)
  bulm_beta_out$model = "bulm_beta"     

  bulm_uncond_out <- get_stan_summaries(y=y, X=X, Psi=Psi, weights=scale_weights_uncond,
                                        n_chains=2, mcmc_burn=mcmc_burn,
                                        mcmc_iter=mcmc_iter, alpha=alpha,
                                        grouped_pop_df=acs_pop_grouped, 
                                        X_formula=X_formula, Psi_formula=Psi_formula)
  bulm_uncond_out$model = "bulm_uncond"     

  # 9. Combine results

  # 8. Fit NPS-informed prior model (Ethan)

  # 9. Fit MRP (Qianyu)



  mrp <- getMRP(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop
  )
  #  mrp_r
  mrpr <- mrp$puma_summary_mrpr

  #  mrp_p
  mrpp <- mrp$puma_summary_mrpp


  mrp1 <- getMRP_INT(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod = modINT
  )


  #  mrp_int
  mrpint <- mrp1$puma_summary_mrpp




  # 10. Fit VSW method (Qi)
  result_VSW <- vsw_out(ps, nps, X_formula) # a vector of 4, (mse, mab, cr, is)
  # 11. Combine results

  results[[sim]] <- rbind(
    direst,
    bulm_ps_out,
    bulm_beta_out,
    bulm_uncond_out,
    result_VSW,
    mrpr,
    mrprp
  )
}

# Aggregate across simulations
results_df <- results %>% list_rbind(names_to = "sim_num")

# Summaries for each method
summary_df <- true_values %>%
  left_join(results_df, by = "PUMA") %>%
  group_by(model) %>%
  summarize(
    MSE = mean((response_true - point_est)^2),
    MAB = mean(abs(response_true - point_est)),
    Coverage = mean(between(response_true, lower_CI, upper_CI)),
    `Int. Score` = mean(int_score(alpha, response_true, lower_CI, upper_CI))
  )

# VSW method doesn't have uncertainty quantification, so I return NA's for them. That'swhy I calculated it alone.

# Save results
save(
  list = c("prob_samples", "nonprob_samples", "results_df", "summary_df"),
  file = "data/ACS_NPS_simulation_results.RData"
)
