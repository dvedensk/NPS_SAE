get_PS <- function(pop.df, type="PPS", samp.frac=.01){#type is PPS or some other size variable
    sample.size <- floor(nrow(pop.df)*samp.frac)
    #create the size variable that will be sampled in proportion to 
    #(making this variable a function of the survey weights induces an informative design)
    size.var <- as.numeric(exp(.2*scale(pop.df$WAGP) + .4*scale(pop.df$PWGTP)))
    inclusion_probs <- inclusionprobabilities(size.var, sample.size)
    inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample.size
    #survey weights are inverse probabilities of selection
    weights <- 1/inclusion_probs
    #Draw Poisson sample: https://en.wikipedia.org/wiki/Poisson_sampling
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
    weights <- weights[sample.idx]

    #for assessing data defect index, we may want to calculate the following: 
    #\rho_{R,G}:
    #cor(pop.df[sample.idx,]$WAGP, popWeights[sample.idx])
    #problem difficulty \sigma2_G:
    #var(pop.df[sample.idx,]$WAGP)

    #return weights since we may use them for calculating properties of the sample, but
    #we will not use them as...
    return(list(weights=weights, idx=sample.idx))
}
