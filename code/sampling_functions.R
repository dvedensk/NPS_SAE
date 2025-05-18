get_PS <- function(pop_df, type="PPS", samp_frac=.01){#type is PPS or some other size variable
    sample_size <- floor(nrow(pop.df) * samp_frac)
    #create the size variable that will be sampled in proportion to 
    #(making this variable a function of the survey weights induces an informative design)
    size_var <- as.numeric(exp(.2*scale(pop.df$WAGP) + .4*scale(pop.df$PWGTP)))
    inclusion_probs <- inclusionprobabilities(size_var, sample_size)
    inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample_size
    #survey weights are inverse probabilities of selection
    weights <- 1/inclusion_probs
    #Draw Poisson sample: https://en.wikipedia.org/wiki/Poisson_sampling
    sample_idx <- which(UPpoisson(inclusion_probs) == 1)
    sample_size <- length(sample_idx)
    weights <- weights[sample_idx]

    return(list(weights=weights, idx=sample_idx))
}

get_NPS <- function(pop_df, noise_level=2, #sd of noise for governing ddi
                    samp_frac=.2, include_internet=F){
    if(!include_internet){
        pop_df <- filter(pop_df, ACCESSINET < 3) #3 indicates no access
    }
    sample_size <- floor(nrow(pop_df)*samp_frac)
    #if R \propto G, ddi = 1, so start there and add noise
    size_var <- as.numeric(1000 + scale(-pop.df$WAGP)) + rnorm(n=nrow(pop_df),
                                                               mean=0,
                                                               sd=noise_level)
    inclusion_probs <- inclusionprobabilities(size_var, sample_size)
    inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample_size
    weights <- 1/inclusion_probs
    sample_idx <- which(UPpoisson(inclusion_probs) == 1)
    sample_size <- length(sample_idx)
    weights <- weights[sample_idx]

    #for assessing data defect index, we may want to calculate the following: 
    #\rho_{R,G}:
    #cor(pop.df[sample.idx,]$WAGP, popWeights[sample.idx])
    #problem difficulty \sigma2_G:
    #var(pop.df[sample.idx,]$WAGP)

    #return weights since we may use them for calculating properties of the sample, but
    #we will not use them as...
    return(list(weights=weights, idx=sample_idx))
}
