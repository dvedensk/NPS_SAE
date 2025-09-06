# library(rstanarm)
# library(brms)

getMRP=function(MR,Ps,P){
  
  
  # ---- 1)model
  nps_ca=MR
  
  # cmdstanr::install_cmdstan()
  fit3 <- rstanarm::stan_glmer(
    HICOV ~ 1 + AGEP + SEX + RAC1P+ (1|PUMA), 
    data = nps_ca, 
    family = binomial, 
    chains = 2
    # backend = "cmdstanr",
    # threads = threading(4)
    # iter = 2000
  )
  message("model done")
  
  
  # ---- 2) poststratification
  #----  Population   MRP-P:using population(P without PWGTP) , MRP-R using probability sample(P with PWGTP)
  #------  
  # MRP-R 
  #------ 
  
  
  # Ensure types are consistent with training data
  {pop <- Ps %>%
    mutate(
      SEX   = factor(SEX,   levels = levels(nps_ca$SEX)),# align factor levels to the TRAINING data
      AGEP  = factor(AGEP,  levels = levels(nps_ca$AGEP)),
      RAC1P = factor(RAC1P, levels = levels(nps_ca$RAC1P)),
      PUMA  = factor(PUMA,  levels = levels(nps_ca$PUMA))
    )
  
  
  pred_df <- pop %>%
    select(SEX, AGEP, RAC1P, PUMA, any_of("PWGTP"))
  
  keep <- complete.cases(pred_df[, c("AGEP","PUMA","SEX","RAC1P")])
  pred_df <- pred_df[keep, , drop = FALSE]
  
  draws <- 1000
  ep <- posterior_epred(fit3, newdata = pred_df[,  c("AGEP","PUMA","SEX","RAC1P")],
                        draws = draws)   # matrix [draws x N_people]
  
  
  # 1 Prep grouping and weights
  puma_vec <- droplevels(factor(pred_df$PUMA)) #droplevels removes any unused levels
  w <- if ("PWGTP" %in% names(pred_df)) pred_df$PWGTP else rep(1, nrow(pred_df))
  w[!is.finite(w) | w < 0] <- 0
  
  
  idx_split <- split(seq_len(nrow(pred_df)), puma_vec)  # list of row indices per PUMA
  
  # For each PUMA, compute weighted mean of ep across rows in that PUMA, for each draw
  puma_draws <- lapply(idx_split, function(idx) {
    wi <- w[idx]
    sw <- sum(wi)
    if (sw == 0) {
      rep(NA_real_, nrow(ep))
    } else {
      wi <- wi / sw
      # ep[, idx] %*% wi # gives a vector of length 'draws'
      as.numeric(ep[, idx, drop = FALSE] %*% wi)
    }
  })
  
  
  # posterior samples and summaries
  puma_draws <- do.call(cbind, puma_draws)
  colnames(puma_draws) <- names(idx_split)
  
  puma_summary_mrpr <- apply(puma_draws, 2, function(d) {
    qs <- quantile(d, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    
    c(point_est = mean(d, na.rm = TRUE),
      sd   = sd(d,   na.rm = TRUE),
      lower_CI  = qs[1],
      upper_CI  = qs[2])
  }) |> t() |> as.data.frame()
  names(puma_summary_mrpr)[names(puma_summary_mrpr) == "lower_CI.2.5%"]  <- "lower_CI"
  names(puma_summary_mrpr)[names(puma_summary_mrpr) == "upper_CI.97.5%"] <- "upper_CI"
  
  
  puma_summary_mrpr$PUMA <- colnames(puma_draws)
  puma_summary_mrpr$model <-"mrp-r"
  }
  
  #------  
  # MRP-p 
  #------
  {
    # Ensure types are consistent with training data
    pop <- P %>%
      mutate(
        SEX   = factor(SEX,   levels = levels(nps_ca$SEX)),# align factor levels to the TRAINING data
        AGEP  = factor(AGEP,  levels = levels(nps_ca$AGEP)),
        RAC1P = factor(RAC1P, levels = levels(nps_ca$RAC1P)),
        PUMA  = factor(PUMA,  levels = levels(nps_ca$PUMA))
      )
    pop$PWGTP=NULL
    
    pred_df <- pop %>%
      select(SEX, AGEP, RAC1P, PUMA, any_of("PWGTP"))
    
    keep <- complete.cases(pred_df[, c("AGEP","PUMA","SEX","RAC1P")])
    pred_df <- pred_df[keep, , drop = FALSE]
    
    draws <- 100
    ep <- posterior_epred(fit3, newdata = pred_df[,  c("AGEP","PUMA","SEX","RAC1P")],
                          draws = draws)   # matrix [draws x N_people]
    
    
    # Prep grouping and weights
    puma_vec <- droplevels(factor(pred_df$PUMA)) #droplevels removes any unused levels
    w <- if ("PWGTP" %in% names(pred_df)) pred_df$PWGTP else rep(1, nrow(pred_df))
    w[!is.finite(w) | w < 0] <- 0
    
    
    idx_split <- split(seq_len(nrow(pred_df)), puma_vec)  # list of row indices per PUMA
    
    puma_draws <- lapply(idx_split, function(idx) {
      wi <- w[idx]
      sw <- sum(wi)
      if (sw == 0) {
        rep(NA_real_, nrow(ep))
      } else {
        wi <- wi / sw
        as.numeric(ep[, idx, drop = FALSE] %*% wi)
      }
    })
    puma_draws <- do.call(cbind, puma_draws)
    colnames(puma_draws) <- names(idx_split)
    
    puma_summary_mrpp <- apply(puma_draws, 2, function(d) {
      qs <- quantile(d, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      c(point_est = mean(d, na.rm = TRUE),
        sd   = sd(d,   na.rm = TRUE),
        lower_CI  = qs[1],
        upper_CI  = qs[2])
    }) |> t() |> as.data.frame()
    
    names(puma_summary_mrpp)[names(puma_summary_mrpp) == "lower_CI.2.5%"]  <- "lower_CI"
    names(puma_summary_mrpp)[names(puma_summary_mrpp) == "upper_CI.97.5%"] <- "upper_CI"
    
    
    puma_summary_mrpp$PUMA <- colnames(puma_draws)
    puma_summary_mrpp$model <-"mrp-p"
    
  }
  
  
  return(
    list(puma_summary_mrpp=puma_summary_mrpp,
         puma_summary_mrpr=puma_summary_mrpr))
  
  
  
  
}
