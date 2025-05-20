require(readr)
require(dplyr)
require(tidyverse)
require(sampling)
require(mvtnorm)
require(survey)
require(purrr)

source(file.path("code","sampling_functions.R"))
source(file.path("code","utils.R"))
source(file.path("code","models", "bulm.R"))

set.seed(99)

#load population file
acs_pop <- read_csv(file.path("data","ACS_NPS_pop.csv")) %>%
    mutate(AGEP = factor(AGEP),
           RAC1P = factor(RAC1P),
           SEX = factor(SEX),
           PUMA = factor(PUMA))

#population values to compare against
true_values <- acs_pop %>% group_by(PUMA) %>%
                           summarize(HICOV=mean(HICOV), 
                                     WAGP=median(WAGP))

acs_pop_grouped <- acs_pop %>%
    group_by(PUMA, AGEP, RAC1P, SEX) %>%
    tally
    
X_formula <- as.formula("~ AGEP + RAC1P + SEX")
Psi_formula <- as.formula("~ -1 + PUMA")
X_pop <-  model.matrix(X_formula , data=acs_pop)
Psi_pop <- model.matrix(Psi_formula, data=acs_pop)

alpha <- .05

##take samples
Nsim <- 100
prob_samples <- nonprob_samples <- results <- list()
for(sim in 1:Nsim) {
  print(sim)
  #do we want to be able to make PS and NPS disjoint?
  ps <- prob_samples[[sim]] <- get_strat_PS(pop_df=acs_pop, samp_frac=.002)
  nps <- nonprob_samples[[sim]] <- get_NPS(pop_df=acs_pop, noise_level=2,
                                           samp_frac=.1, include_internet=F)

  ps_scale_weights <- length(ps$idx) *  ps$weights/sum(ps$weights)
  nps_scale_weights <- length(nps$idx) *  nps$weights/sum(nps$weights)

  ps <- acs_pop[ps$idx, ]
  nps <- acs_pop[nps$idx, ]

  X_ps <- model.matrix(X_formula , data=ps)
  Psi_ps <- model.matrix(Psi_formula, data=ps)
  y_ps <- ps$HICOV

  # Fit models here.
  # Each model should output a data frame with one row for each PUMA
  # The columns should be (1) PUMA ID (2) point estimate, (3) lower and (4) upper
  # endpoints of a CI with confidence level alpha, (4) the method name
  samp.design <- svydesign(ids=~1, weights=~PWGTP, data=ps)
  direst <- svyby(~HICOV, ~ PUMA, samp.design, svymean, vartype="se")               
  direst <- direst %>% rename(point_est = HICOV) %>%
                       mutate(lower_CI = point_est + qnorm(alpha/2)*se,
                              upper_CI = point_est + qnorm(1-alpha/2)*se) %>%
                       select(-se)
  rownames(direst) <- c()
  direst$model <- "direst"
       
  bulm_out <- bulm_results(grouped_pop_df = acs_pop_grouped,
                           X = X_ps,
                           Psi = Psi_ps,
                           y = y_ps,
                           weights = ps_scale_weights,
                           sigma2_beta = 1e4,
                           iter = 2000,
                           burn = 1000,
                           alpha = alpha,
                           X_formula = X_formula,
                           Psi_formula = Psi_formula,
                           summaries=TRUE)
  bulm_out$model <- "bulm"

  # combine models into a single data.frame here
  results[[sim]] <-rbind(direst,
                         bulm_out)
}
results_df <- results %>% list_rbind(names_to="sim_num")

summary_df <- true_values %>%
    left_join(results_df, by="PUMA") %>%
    group_by(model) %>%
    summarize(MSE = mean( (HICOV - mean)^2 ),
              MAB = mean(abs(HICOV - mean(mean))),
              Coverage = mean(between(HICOV, lower_CI, upper_CI)),
              `Int. Score` = mean(int_score(alpha, HICOV, lower_CI, upper_CI)))
              
save(list(prob_samples=prob_samples,
          nonprob_samples=nonprob_samples,
          results_df=results_df,
          summary_df=summary_df),
     file="data/ACS_NPS_simulation_results.RData")
