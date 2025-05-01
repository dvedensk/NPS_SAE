require(readr)
require(dplyr)
require(tidyverse)
require(sampling)
require(mvtnorm)

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "nps_prior.R"))

set.seed(99)

#load population file
read_csv(file.path("data", "ACS_NPS_pop.csv"))

#population values to compare against
true.values <- acs.pop %>% group_by(PUMA) %>%
                           summarize(HICOV=mean(HICOV), 
                                     WAGP=median(WAGP))

##take samples
Nsim <- 100
prob.samples <- non.prob.samples <- list()
for(sim in 1:Nsim){
    #do we want to be able to make PS and NPS disjoint?
    ps <- prob.samples[[sim]] <- get_PS(pop.df=acs.pop, samp.frac=.002)
    nps <- non.prob.samples[[sim]] <- get_NPS(pop.df=acs.pop, noise.level=2,
                                              samp.frac=.1, include.internet=F)

    ps.scale.weights <- length(ps$idx) *  ps$weights/sum(ps$weights)
    nps.scale.weights <- length(nps$idx) *  nps$weights/sum(nps$weights)

    ps <- acs.pop[ps$idx, ]
    nps <- acs.pop[nps$idx, ]
    #fit models here:

}

save(list(prob.samples=prob.samples,
          non.prob.samples=non.prob.samples),
     file=file.path("data", "ACS_NPS_samples.RData"))
     