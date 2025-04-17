library(readr)
library(dplyr)
# library(tidyverse) # faster, safer, and requires less memory to 
                     # import individual packages from the tidyverse
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
households <- read_csv("psam_h06.csv", col_select=c("SERIALNO", "ACCESSINET", "TEL"))
persons <- read_csv("psam_p06.csv", col_select=c("SERIALNO",
                                                 "PWGTP",
                                                 "PUMA",
                                                 "STATE",
                                                 "AGEP",
                                                 "RAC1P", #Race (9 levels)
                                                 "SEX",
                                                 "WAGP", #Wages/salary past 12 months
                                                 "SCHL", #Collapse into fewer categories
                                                 "HICOV")) #Health insurance (1 = covered, 2 = not covered)

acs.pop <- persons %>% left_join(households, by=c("SERIALNO"))
acs.pop <- acs.pop %>% filter(substr(SERIALNO, start=5,stop=6)=="HU") %>% #exclude group quarters
                       filter(!is.na(WAGP))

write_csv(acs.pop, file="ACS_NPS_pop.csv")

##take samples
Nsim <- 100
prob.samples <- non.prob.samples <- list()
for(sim in 1:Nsim){
    #do we want to be able to make PS and NPS disjoint?
    prob.samples[[sim]] <- ps <- get_PS(pop.df=acs.pop, samp.frac=.002)
    non.prob.samples[[sim]] <- nps <- get_NPS(pop.df=acs.pop, noise.level=2,
                                              samp.frac=.1, include.internet=F)

    ps.scale.weights <- length(ps$idx) *  ps$weights/sum(ps$weights)
    nps.scale.weights <- length(nps$idx) *  nps$weights/sum(nps$weights)

    ps <- acs.pop[ps$idx, ]
    nps <- acs.pop[nps$idx, ]
    #fit models here:
}

all.samples <- list(prob.samples=prob.samples,
                    non.prob.samples=non.prob.samples)

save(all.samples,
     file="ACS_NPS_samples.RData")

# TODO: test on small data w/ noninformative prior
nps_post <- function(y, X, mu_beta, k0, V, a, b) {
    n <- nrow(X)
    p <- ncol(X)

    XTX <- crossprod(X)
    inv_scale_mat <- solve(XTX + solve(k0 * V)) # TODO: check why eqn. (9) is different
    
    beta_hat <- solve(XTX) %*% crossprod(X, y) # eqn 7
    W <- inv_scale_mat %*% XTX # eqn 6
    post_mn <- W %*% beta_hat + (diag(1, p) - W) %*% mu_beta # eqn 5

    RSS <- crossprod(y - X %*% beta_hat) # eqn 10
    bet_minus_mu <- beta_hat - mu_beta
    SS <- RSS + crossprod(bet_minus_mu, solve(solve(XTX) + k0 * V) %*% bet_minus_mu) # eqn 9
    post_var <- inv_scale_mat %*% (SS + 2 * b) / (n + 2 * a - 2) # eqn 8

    return(list(mn = post_mn, vr = post_var))
}
