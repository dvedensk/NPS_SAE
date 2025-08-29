vsw_out <- function(ps, nps, X_formula){
  # 3. Build design matrices for PS
  X_ps   <- model.matrix(X_formula,   data = ps)
  y_ps   <- ps$HICOV
  X_nps <- model.matrix(X_formula, data = nps)
  y_nps <- nps$HICOV
  
  
  
  #fit models here:
  model_ps <- glm(y_ps~X_ps, family = binomial(link = "logit"))
  model_nps <- glm(y_nps~X_nps, family = binomial(link = "logit"), data = nps)
  Z_hat_ps <- data.frame(
    "PUMA" = ps$PUMA, "prediction" = predict(model_ps, type = "response")
  ) %>% group_by(PUMA) %>% summarise(Z_k1 = mean(prediction), Z_k0 = 1-mean(prediction)) 
  
  Z_hat_nps <- data.frame(
    "PUMA" = nps$PUMA, "prediction" = predict(model_nps, type = "response")
  ) %>% group_by(PUMA) %>% summarise(Z_k1 = mean(prediction), Z_k0 = 1-mean(prediction)) 
  
  Z_combined <- Z_hat_ps %>%
    left_join(Z_hat_nps, by = "PUMA", suffix = c("_ps", "_nps"))
  
  nps_count <- nps %>%
    count(PUMA, name = "n_nps")%>%
    semi_join(ps, by = "PUMA")
  
  ps_count <- ps %>% count(PUMA, name = "n_nps")
  
  K <- nrow(Z_hat_ps)
  C <- length(unique(ps$HICOV))
  # if(nrow(Z_hat_ps) != nrow(Z_hat_nps)){stop("Probability sample and non-probability sample domains dismatch.")}
  hat_beta_c <- colMeans(Z_combined[,4:5] - Z_hat_ps[,2:3])
  hat_sig <- sum((Z_combined[,4:5] - Z_combined[,2:3]) - 
                   matrix(hat_beta_c, nrow = nrow(Z_combined), ncol = 2, byrow = TRUE)^2)/(C*(K-1))
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
  
  mse <- mean((combined_results[,1] - true_values$HICOV)^2)
  mab <- mean(abs(combined_results[,1] - true_values$HICOV))
  cr = is <- rep(NA, length(mse))
  summary_df <- cbind(mse, mab, cr, is)
  
  return(summary_df)
}
  




