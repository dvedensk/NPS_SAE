# library(readr) # Redundant with tidyverse
# library(dplyr) # Redundant with tidyverse
library(tidyverse)
library(sampling)
library(mvtnorm)
library(survey)
# library(purrr) # Redundant with tidyverse
library(Matrix)
library(LaplacesDemon)
library(parallelly)
library(cmdstanr)
library(rmarkdown)

models_path <- file.path("code", "models")

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R"))
source(file.path(models_path, "bulm.R"))
source(file.path(models_path, "VSW.R"))

# Render the Rmd file in the global environment 
# (effectively sources the R code chunks)
render(
  file.path(models_path, "nps_prior_logistic.Rmd"), 
  output_format = "html_document", 
  output_file = "nps_prior_logistic.html",
  output_dir = models_path,
  envir = globalenv(), 
  quiet = TRUE
) %>% 
  suppressWarnings() # Don't throw deprecation warnings

source(file.path(models_path, "mrp_all.R"))
# load .stan
mod <- cmdstan_model(
  file.path(models_path, "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

modINT <- cmdstan_model(
  file.path(models_path, "mrp_int2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

Sys.setenv(STAN_NUM_THREADS = availableCores() - 1) # number of threads should always be, at most, one fewer than the number of available cores (leave one for system processes.) 
# Also, parallelly::availableCores() - 1 is safer (for example, on HPC) because it also fulfills SLURM constraints. 
# If we do this on HPC, we can safely use availableCores() instead of availableCores() - 1.

set.seed(99)

# ============================================================
# CONFIGURATION: Response variable and sampling weights
# ============================================================
# Response variable to analyze (must exist in ACS_NPS_pop.csv)
response_var <- "PUBCOV" # Default: PUBCOV (binary). Can also use HICOV, WAGP, etc.

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
  ps <- get_strat_PS(pop_df = acs_pop, samp_frac = .005, weight_config = PS_weight_config)
  nps <- get_NPS(pop_df = acs_pop, samp_frac = .05, weight_config = NPS_weight_config, internet_only = FALSE)

  prob_samples[[sim]] <- ps
  nonprob_samples[[sim]] <- nps

  # Calculate and print sample diagnostics (DDC and ESS)
  diagnostics <- calculate_sample_diagnostics(
    pop_df = acs_pop,
    ps_sample = ps,
    nps_sample = nps,
    response_var = response_var,
    print_results = TRUE,
    sim_num = sim
  )

  # 2. Scale weights for pseudolikelihood models
  ps_scale_weights <- length(ps$idx) * ps$weights / sum(ps$weights)
  nps_scale_weights <- length(nps$idx) * nps$weights / sum(nps$weights)

  # 3. Extract PS and NPS data frames
  ps <- acs_pop[ps$idx, ]
  nps <- acs_pop[nps$idx, ]

  PUMA <- ps$PUMA
  PUMA_levels <- unique(PUMA)
  
  # 4. Build design matrices for PS and NPS separately
  X_ps <- model.matrix(X_formula, data = ps)
  Psi_ps <- model.matrix(Psi_formula, data = ps)
  y_ps <- ps[[response_var]]

  X_nps <- model.matrix(X_formula, data = nps)
  Psi_nps <- model.matrix(Psi_formula, data = nps)
  y_nps <- nps[[response_var]]

  # 5. Direct estimate on probability sample
  samp.design <- svydesign(ids = ~1, weights = ~PWGTP, data = ps)

  # svyby extraction: handle both binary and continuous responses
  direst_raw <- svyby(
    as.formula(paste0("~", response_var)),
    ~PUMA,
    samp.design,
    svymean,
    vartype = "var",
    na.rm = TRUE
  ) %>% arrange(PUMA)

  # Extract correct columns based on response type (binary vs continuous)
  # Check if response variable is binary (0/1) by examining unique values
  unique_vals <- unique(acs_pop[[response_var]])
  unique_vals <- unique_vals[!is.na(unique_vals)]
  is_binary <- all(unique_vals %in% c(0, 1)) && length(unique_vals) <= 2

  all_names <- names(direst_raw)
  response_cols <- setdiff(all_names, c("PUMA", all_names[grepl("^var", all_names)]))

  if (is_binary && length(response_cols) > 1) {
    # Binary variable with factor levels: select column ending in 1 or TRUE
    response_col <- response_cols[grepl("(1|TRUE)$", response_cols)]
    if (length(response_col) == 0) response_col <- response_cols[1]
    var_col <- paste0("var.", response_col)
  } else {
    # Continuous or numeric binary (0/1): single column
    response_col <- response_cols[1]
    var_col <- "var"
  }

  direst <- direst_raw %>%
    mutate(
      point_est = .data[[response_col]],
      se = sqrt(.data[[var_col]]),
      lower_CI = point_est + qnorm(alpha / 2) * se,
      upper_CI = point_est + qnorm(1 - alpha / 2) * se,
      model = "direst"
    ) %>%
    select(PUMA, point_est, lower_CI, upper_CI, model)



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


  # 7. Estimate IP weights for NPS and fit unit-level model
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

  # 10. Fit NPS-informed prior model (Ethan)

  # TODO: use survey::rake to make the survey weight covariate
  # TODO: check mixing; try other PG sampler packages; if nothing else, use Stan for power prior 

  nps_prior_res <- nps_prior_mcmc(
    y_ps, 
    cbind(X_ps, Psi_ps), 
    y_nps, 
    cbind(X_nps, Psi_nps), 
    niter = 4000, 
    PUMA_levels = PUMA_levels,
    wts = ps_scale_weights, 
    typeIerr = alpha
  )
  
  # 11. Fit MRP (Qianyu)
  mrp=getMRP(MR=nps,
             ps=ps,
             acs_pop=acs_pop)
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
  mrpp=mrp$puma_summary_mrpp[,c("PUMA","point_est","lower_CI","upper_CI","model")]
  
  # 12. Fit VSW method (Qi)
  result_VSW <- vsw_out(ps, nps, X_formula) # a vector of 4, (mse, mab, cr, is)

  # 13. Combine results
  results[[sim]] <- rbind(
    direst,
    bulm_out,
    bulm_ipw,
    result_VSW,
    mrpr,
    mrprp,
    nps_prior_res
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
