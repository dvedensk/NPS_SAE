require(readr)
require(dplyr)
require(tidyverse)
require(sampling)
source("sampling_functions.R")

set.seed(99)

#load population file
read_csv("data/ACS_NPS_pop.csv")

#population values to compare against
true.values <- acs_pop %>% group_by(PUMA) %>%
                           summarize(HICOV=mean(HICOV), 
                                     WAGP=median(WAGP))

##take samples
Nsim <- 100
prob_samples <- non_prob_samples <- list()
for(sim in 1:Nsim){
    #do we want to be able to make PS and NPS disjoint?
    ps <- prob_samples[[sim]] <- get_PS(pop.df=acs_pop, samp_frac=.002)
    nps <- nonprob_samples[[sim]] <- get_NPS(pop_df=acs_pop, noise_level=2,
                                              samp_frac=.1, include_internet=F)

    ps.scale.weights <- length(ps$idx) *  ps$weights/sum(ps$weights)
    nps.scale.weights <- length(nps$idx) *  nps$weights/sum(nps$weights)

    ps <- acs_pop[ps$idx, ]
    nps <- acs_pop[nps$idx, ]
    #fit models here:

}

save(list(prob_samples=prob_samples,
          nonprob.samples=non_prob.samples),
     file="data/ACS_NPS_samples.RData")
