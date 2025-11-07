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

source(file.path("code","sampling_functions.R"))
source(file.path("code","utils.R"))   # utils.R must define estimate_ipw()
source(file.path("code","models","bulm.R"))
source(file.path("code","models","VSW.R"))

set.seed(99)

# load population file
acs_pop <- read_csv(file.path("data","ACS_NPS_pop.csv")) %>%
  mutate(AGEP = factor(AGEP),
         RAC1P = factor(RAC1P),
         SEX = factor(SEX),
         PUMA = factor(PUMA))

# population values to compare against
true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(HICOV = mean(HICOV),
            WAGP  = median(WAGP))

acs_pop_grouped <- acs_pop %>%
  group_by(PUMA, AGEP, RAC1P, SEX) %>%
  tally

X_formula   <- as.formula("~ AGEP + RAC1P + SEX")
Psi_formula <- as.formula("~ -1 + PUMA")

alpha <- .05

# take samples
Nsim <- 100
prob_samples    <- list()
nonprob_samples <- list()
results         <- list()
summary_df_VSW <- list()

sim <- 1

ps  <- get_strat_PS(pop_df = acs_pop, samp_frac = .002, w_wage=.2, w_pwgt=.2, w_cit=.5)
  
ps_df  <- acs_pop[ps$idx, ]
ps_df$weights <- ps$weights

ps_df %>%
    group_by(PUMA) %>%
    summarize(cov=mean(HICOV)) %>%
    ggplot() +
    geom_histogram(aes(x=cov))

plot(table(ps_df$RAC1P, ps_df$SEX))
plot(table(ps_df$RAC1P, ps_df$AGEP))

# Try other filtering variables. Maybe HISPEED?
nps <- get_NPS(pop_df = acs_pop, noise_level = 2,
               samp_frac = .1, include_internet = FALSE)

nps_df  <- acs_pop[nps$idx, ]
nps_df$weights <- nps$weights

#Check design effect

ps_df$source="ps"
nps_df$source="nps"
acs_df <- acs_pop
acs_df$weights = 1
acs_df$source="pop"

rbind(ps_df, nps_df,acs_df) %>% 
    group_by(PUMA, source, RAC1P, SEX) %>%
    summarize(cov=mean(HICOV)) %>%
    ggplot() +
    geom_histogram(aes(x=cov, fill=source), color="black", alpha=.7) +
    facet_wrap(SEX~RAC1P, scales="free")


rbind(ps_df, nps_df) %>% 
    group_by(PUMA, source, AGEP, SEX) %>%
    summarize(cov=mean(HICOV)) %>%
    ggplot() +
    geom_histogram(aes(x=cov, fill=source), color="black", alpha=.7) +
    facet_wrap(SEX~AGEP, scales="free") +
    scale_color_colorblind() +
    theme_bw()

