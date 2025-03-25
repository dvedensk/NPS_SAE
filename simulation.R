library(readr)
library(dplyr)
library(tidyverse)
library(lobstr)
library(sampling)

set.seed(99)

get_PS <- function(pop.df, type="PPS", samp.frac=.01){#type is PPS or some other size variable
    sample.size <- floor(nrow(pop.df)*samp.frac)
    size.var <- as.numeric(exp(.2*scale(pop.df$WAGP) + .4*scale(pop.df$PWGTP)))
    inclusion_probs <- inclusionprobabilities(size.var, sample.size)
    inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample.size
    weights <- 1/inclusion_probs
    sample.idx <- which(UPpoisson(inclusion_probs) == 1)
    sample.size <- length(sample.idx)
    weights <- weights[sample.idx]

    return(list(weights=weights, idx=sample.idx))
} 

get_NPS <- function(pop.df, noise.level=2, #sd of noise for governing ddi
                    samp.frac=.2, include.internet=F){
    if(!include.internet){
        pop.df <- filter(pop.df, ACCESSINET < 3) #3 indicates no access
    }
    sample.size <- floor(nrow(pop.df)*samp.frac)
    #if R \propto G, ddi = 1, so start there and add noise
    size.var <- as.numeric(1000 + scale(-pop.df$WAGP)) + rnorm(n=nrow(pop.df),
                                                               mean=0,
                                                               sd=noise.level)
    inclusion_probs <- inclusionprobabilities(size.var, sample.size)
    inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample.size
    weights <- 1/inclusion_probs
    sample.idx <- which(UPpoisson(inclusion_probs) == 1)
    sample.size <- length(sample.idx)
    
    #for assessing data defect index, we may want to calculate the following: 
    #\rho_{R,G}:
    #cor(pop.df[sample.idx,]$WAGP, popWeights[sample.idx])
    #problem difficulty \sigma2_G:
    #var(pop.df[sample.idx,]$WAGP)

    return(list(weights=weights, idx=sample.idx))    
}

##Process the population file
#filenames <- paste0("psam_pus",c("a","b","c","d"),".csv")
households <- read_csv("psam_husa.csv", col_select=c("SERIALNO", "ACCESSINET", "TEL"))
persons <- read_csv("psam_pusa.csv", col_select=c("SERIALNO",
                                                 "PWGTP",
                                                 "PUMA", 
                                                 "STATE", 
                                                 "AGEP",
                                                 "RAC1P", #Race (9 levels)
                                                 "SEX",
                                                 "WAGP", #Wages/salary past 12 months
                                                 "SCHL", #Collapse into fewer categories
                                                 "HICOV"))

acs.pop <- persons %>% left_join(households, by=c("SERIALNO"))
acs.pop <- acs.pop %>% filter(substr(SERIALNO, start=5,stop=6)=="HU") %>% #exclude group quarters
                       filter(!is.na(WAGP))

write_csv(acs.pop, file="ACS_NPS_pop.csv")

##take samples
Nsim <- 100
prob.samples <- non.prob.samples <- list()
for(sim in 1:Nsim){
    #do we want to be able to make PS and NPS disjoint?
    ps <- prob.samples[[sim]] <- get_PS(pop.df=acs.pop, samp.frac=.01)
    nps <- non.prob.samples[[sim]] <- get_NPS(pop.df=acs.pop, noise.level=2,
                                              samp.frac=.1, include.internet=F)
    
    ps.scale.weights <- length(ps$idx) *  ps$weights/sum(ps$weights)
    nps.scale.weights <- length(nps$idx) *  nps$weights/sum(nps$weights)

    #fit models here:
}

save(list(prob.samples=prob.samples,
          non.prob.samples=non.prob.samples),
     file="ACS_NPS_samples.RData"))
