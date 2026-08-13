vsw_out <- function(ps, nps, X_formula = NULL, response){
  # 3. Build design matrices for PS
  # X_ps   <- model.matrix(X_formula,   data = ps)
  # X_nps <- model.matrix(X_formula, data = nps)
  y_ps  <- ps[[response]]
  y_nps <- nps[[response]]
  
  all_dat <- rbind(ps, nps)
  # X_all <- model.matrix(X_formula, data = all_dat)
  y_all <- all_dat[[response]]
  
  
  #fit models here:
  # model_ps <- glm(y_ps~X_ps, family = binomial(link = "logit"))
  # model_nps <- glm(y_nps~X_nps, family = binomial(link = "logit"), data = nps)
  # model_all <- glm(y_all~X_all, family = binomial(link = "logit"), data = all_dat)
  
  z_hat_pool <- data.frame(
    "PUMA" = all_dat$PUMA, "prediction" = all_dat[[response]]
  ) %>% group_by(PUMA) %>% summarise(Z_k1 = mean(prediction), Z_k0 = 1-mean(prediction))
  
  
  Z_hat_ps <- data.frame(
    "PUMA" = ps$PUMA, "prediction" = ps[[response]]
  ) %>% group_by(PUMA) %>% summarise(Z_k1 = mean(prediction), Z_k0 = 1-mean(prediction)) 
  
  Z_hat_nps <- data.frame(
    "PUMA" = nps$PUMA, "prediction" = nps[[response]]
  ) %>% group_by(PUMA) %>% summarise(Z_k1 = mean(prediction), Z_k0 = 1-mean(prediction)) 
  
  Z_combined <- Z_hat_ps %>%
    left_join(Z_hat_nps, by = "PUMA", suffix = c("_ps", "_nps"))
  
  nps_count <- nps %>%
    dplyr::count(PUMA, name = "n_nps")%>%
    semi_join(ps, by = "PUMA")
  
  ps_count <- ps %>% dplyr::count(PUMA, name = "n_ps")
  
  K <- nrow(Z_hat_ps)
  C <- length(unique(ps[[response]]))
  # if(nrow(Z_hat_ps) != nrow(Z_hat_nps)){stop("Probability sample and non-probability sample domains mismatch.")}
  hat_beta_c <- colMeans(Z_combined[,4:5] - Z_hat_ps[,2:3])
  hat_sig <- sum((Z_combined[,4:5] - Z_combined[,2:3] - 
                    matrix(hat_beta_c, nrow = nrow(Z_combined), ncol = 2, byrow = TRUE))^2)/(C*(K-1))
  hat_nu_kc <- Z_combined[,4:5]*(1-Z_combined[,4:5])
  combined_results <- matrix(NA, nrow = K, ncol = C)
  for (k in 1:K) {
    for(c in 1:C){
      curr_emse_np <- as.numeric(hat_sig + hat_nu_kc[k,c]/(nps_count[k,2]-1) + hat_beta_c[c]^2)
      curr_emse_p <- as.numeric(1/ps_count[k,2]*(
        nps_count[k,2]/(nps_count[k,2]-1)*hat_nu_kc[k,c] +
          hat_beta_c[c]*(2*Z_combined[k,c+3]-1) - hat_beta_c[c]^2 - hat_sig
      ))
      curr_w <- curr_emse_np/(curr_emse_np+curr_emse_p)
      combined_results[k,c] <- as.numeric(curr_w*Z_combined[k,c+1] + (1-curr_w)*Z_combined[k,c+3])
    }
  }
  
  return(
    # This version only allow the case when the data has only two categores, and the function is estimating the probability of the response 1.
    data.frame("PUMA" = Z_hat_nps$PUMA, "VSW_point_est" = combined_results[,1], "ps_est" = Z_hat_ps[,2,drop = TRUE], "nps_est" = Z_hat_nps[,2,drop = TRUE], "pooled_results" = z_hat_pool[,2, drop = TRUE], "lower_CI" = NA, "upper_CI" = NA, "model" = "VSW") )
}
