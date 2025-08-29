int_score <- function(alpha, truth, L, U){
  return(
    (U - L) + 2/alpha*(truth < L)*(L - truth) + 2/alpha*(truth > U)*(truth - U)
  )
}

# Estimate pseudo‐inclusion weights for a non‐probability sample
# via a logistic model distinguishing prob vs. non‐prob cases

estimate_ipw <- function(ps, nps, cov_formula) {
  library(dplyr)
  # 1. Stack and label
  combined <- bind_rows(
    ps  %>% mutate(Z = 0),
    nps %>% mutate(Z = 1)
  )
  
  # 2. Fit propensity‐to‐be‐nonprob logistic regression
  fit_prop <- glm(
    formula = as.formula(paste0("Z ", deparse(cov_formula))),
    family  = binomial(),
    data    = combined
  )
  
  # 3. Predict P(Z=1 | X) on nonprob units
  p_hat  <- predict(fit_prop, newdata = nps, type = "response")
  ipw_raw <- 1 / p_hat
  
  # 4. Normalize so sum(ipw) = nrow(nps) (keeps scale comparable to PS size)
  ipw <- ipw_raw * (nrow(nps) / sum(ipw_raw))
  return(ipw)
}
