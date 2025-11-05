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

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R")) # utils.R must define estimate_ipw()
source(file.path("code", "models", "bulm.R"))
source(file.path("code", "models", "VSW.R"))

source(file.path("code", "models", "mrp_all.R"))
# load .stan
mod <- cmdstan_model(
  file.path("code", "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

modINT <- cmdstan_model(
  file.path("code", "mrp_int2.stan"),
  cpp_options = list(stan_threads = TRUE)
)
Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())

set.seed(99)

# load population file
acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv")) %>%
  mutate(
    AGEP = factor(AGEP),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# population values to compare against
true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    HICOV = mean(HICOV),
    WAGP = median(WAGP)
  )

acs_pop_grouped <- acs_pop %>%
  group_by(PUMA, AGEP, RAC1P, SEX) %>%
  tally()

X_formula <- as.formula("~ AGEP + RAC1P + SEX")
Psi_formula <- as.formula("~ -1 + PUMA")

alpha <- .05

# take samples
Nsim <- 1
prob_samples <- list()
nonprob_samples <- list()
results <- list()
summary_df_VSW <- list()
for (sim in 1:Nsim) {
  print(sim)

  # 1. Draw probability and nonprobability samples
  ps <- get_strat_PS(pop_df = acs_pop, samp_frac = .002)
  nps <- get_NPS(
    pop_df = acs_pop, noise_level = 2,
    samp_frac = .1, include_internet = FALSE
  )

  prob_samples[[sim]] <- ps
  nonprob_samples[[sim]] <- nps

  # 2. Scale weights for pseudolikelihood models
  ps_scale_weights <- length(ps$idx) * ps$weights / sum(ps$weights)
  nps_scale_weights <- length(nps$idx) * nps$weights / sum(nps$weights)

  # 3. Extract PS and NPS data frames
  ps <- acs_pop[ps$idx, ]
  nps <- acs_pop[nps$idx, ]

  # 4. Build design matrices for PS and NPS separately
  X_ps <- model.matrix(X_formula, data = ps)
  Psi_ps <- model.matrix(Psi_formula, data = ps)
  y_ps <- ps$HICOV

  X_nps <- model.matrix(X_formula, data = nps)
  Psi_nps <- model.matrix(Psi_formula, data = nps)
  y_nps <- nps$HICOV

  # 5. Direct estimate on probability sample
  samp.design <- svydesign(ids = ~1, weights = ~PWGTP, data = ps)
  direst <- svyby(~HICOV, ~PUMA, samp.design, svymean, vartype = "se") %>%
    rename(point_est = HICOV) %>%
    mutate(
      lower_CI = point_est + qnorm(alpha / 2) * se,
      upper_CI = point_est + qnorm(1 - alpha / 2) * se
    ) %>%
    select(-se) %>%
    mutate(model = "direst")


  # 5. BULM on probability sample
  # Sho: what is the difference between 5 and 6 exactly?

  # 6. Fit unit-level model on probability sample

  bulm_out <- bulm_results(
    grouped_pop_df = acs_pop_grouped,
    alpha          = alpha,
    X              = X_ps,
    Psi            = Psi_ps,
    y              = y_ps,
    weights        = ps_scale_weights,
    sigma2_beta    = 1e4,
    iter           = 2000,
    burn           = 1000,
    X_formula      = X_formula,
    Psi_formula    = Psi_formula,
    summaries      = TRUE
  )
  bulm_out$model <- "bulm"

  # 7. Estimate IP weights for NPS
  ipw <- estimate_ipw(ps = ps, nps = nps, cov_formula = X_formula)

  # 9. BULM on nonprobability sample with IPW


  # 10. Estimate IP weights for NPS and fit unit-level model
  #    on nonprobability sample with them (Tracy)
  ipw <- estimate_ipw(ps = ps, nps = nps, cov_formula = X_formula)


  bulm_ipw <- bulm_results(
    grouped_pop_df = acs_pop_grouped,
    alpha          = alpha,
    X              = X_nps,
    Psi            = Psi_nps,
    y              = y_nps,
    weights        = ipw,
    sigma2_beta    = 1e4,
    iter           = 2000,
    burn           = 1000,
    X_formula      = X_formula,
    Psi_formula    = Psi_formula,
    summaries      = TRUE
  )
  bulm_ipw$model <- "bulm_ipw"

  # 11. Fit NPS-informed prior model (Ethan)

  # 12. Fit MRP (Qianyu)
  # ------ MRP-P and MRP-R ------
  mrp <- getMRP(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop
  )

  mrpr <- mrp$puma_summary_mrpr
  mrpp <- mrp$puma_summary_mrpp

  # ------ MRP-INT ------
  mrp1 <- getMRP_INT(
    MR = nps,
    ps = ps,
    acs_pop = acs_pop,
    mod = modINT
  )
  mrpint <- mrp1$puma_summary_mrpp


  # 13. Fit VSW method (Qi)
  result_VSW <- vsw_out(ps, nps, X_formula) # a vector of 4, (mse, mab, cr, is)

  # 14. Combine results
  results[[sim]] <- rbind(
    direst,
    bulm_out,
    bulm_ipw,
    result_VSW,
    mrpr,
    mrpp,
    mrpint
  )
}

# Aggregate across simulations
results_df <- results %>% list_rbind(names_to = "sim_num")

# Summaries for each method
summary_df <- true_values %>%
  left_join(results_df, by = "PUMA") %>%
  group_by(model) %>%
  summarize(
    MSE = mean((HICOV - point_est)^2),
    MAB = mean(abs(HICOV - point_est)),
    Coverage = mean(between(HICOV, lower_CI, upper_CI)),
    `Int. Score` = mean(int_score(alpha, HICOV, lower_CI, upper_CI))
  )

# VSW method doesn't have uncertainty quantification, so I return NA's for them. That'swhy I calculated it alone.

# Save results
save(
  list = c("prob_samples", "nonprob_samples", "results_df", "summary_df"),
  file = "data/ACS_NPS_simulation_results.RData"
)
