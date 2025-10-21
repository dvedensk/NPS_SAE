




getMRP=function(MR=nps,
                ps=ps,
                acs_pop=acs_pop){
  
  # TRAINING DATA ----
  nps_ca <- MR  # nps
  
  # lock factor levels from training
  nps_ca$AGEP  <- factor(nps_ca$AGEP)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev <- levels(factor(nps_ca$PUMA))  # lock PUMA levels
  
  # fixed-effects design 
  X_train <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train); n <- nrow(X_train)
  y <- as.integer(nps_ca$HICOV)
  group <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  J <- length(PUMA_lev)
  
  # Helper to build prediction matrices that EXACTLY match training columns
  mk_pred <- function(df, X_train_cols, PUMA_lev, train_ref) {
    if (is.null(df) || nrow(df) == 0) {
      # return empty structures with correct shapes/types
      return(list(
        N = 0,
        X = matrix(0, 0, length(X_train_cols), dimnames = list(NULL, X_train_cols)),
        g = integer(0),
        df_out = df[0, , drop = FALSE]
      ))
    }
    df$AGEP  <- factor(df$AGEP,  levels = levels(train_ref$AGEP))
    df$SEX   <- factor(df$SEX,   levels = levels(train_ref$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(train_ref$RAC1P))
    df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
    
    Xp <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = df)
    # add any missing training columns (if some levels absent in this df)
    miss <- setdiff(X_train_cols, colnames(Xp))
    if (length(miss)) {
      Xp <- cbind(Xp, matrix(0, nrow(Xp), length(miss),
                             dimnames = list(NULL, miss)))
    }
    # order columns to match training
    Xp <- Xp[, X_train_cols, drop = FALSE]
    
    list(
      N  = nrow(Xp),
      X  = Xp,
      g  = as.integer(df$PUMA),
      df_out = df
    )
  }
  
  # Build Ps (probability-sample poststrat) block
  ps_b <- mk_pred(ps, colnames(X_train), PUMA_lev, nps_ca)  # ps include PWGTP for your weighting step
  w_ps <- ps$PWGTP[match(rownames(ps_b$X), rownames(ps))]  # align if needed
  
  # Build P (full population poststrat) block
  p_b  <- mk_pred(acs_pop,  colnames(X_train), PUMA_lev, nps_ca)
  w_pop <- rep(1, p_b$N)
  
  
  
  stan_data <- list(
    # training
    n = n,
    k = k,
    X = X_train,
    y = y,
    J = J,
    group = group,
    grainsize = 500      # choose 250–1000; tune for speed
    # 
    # # Ps (probability-sample poststratification)
    # Nps = ps_b$N,
    # Xps = ps_b$X,
    # group_ps = ps_b$g,
    # 
    # # P (full-population poststratification)
    # Npop = p_b$N,
    # Xpop = p_b$X,
    # group_pop = p_b$g
  )
  
  
  stan_data <- c(stan_data, list(
    Nps = ps_b$N, Xps = ps_b$X, group_ps = ps_b$g, w_ps = as.vector(w_ps),
    Npop = p_b$N, Xpop = p_b$X, group_pop = p_b$g, w_pop = as.vector(w_pop)
  ))
  
  
  fit <- mod$sample(
    data = stan_data,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  ) #176.1
  
  
  
  p_ps_draws  <- fit$draws("mu_ps",  format = "matrix")  # S x Nps
  p_pop_draws <- fit$draws("mu_pop", format = "matrix")  # S x Npop
  
  
  
  
  make_summary <- function(draw_mat, puma_names, tag) {
    draw_mat[!is.finite(draw_mat)] <- NA_real_
    tibble(PUMA = puma_names) |>
      bind_cols(as.data.frame(t(apply(draw_mat, 2, function(d) {
        qs <- quantile(d, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
        c(point_est = mean(d, na.rm = TRUE),
          sd        = sd(d,   na.rm = TRUE),
          lower_CI  = qs[1],
          upper_CI  = qs[2])
      })))) |>
      mutate(model = tag, .before = 1)
  }
  
  
  
  mu_ps_draws  <- fit$draws("mu_ps",  format = "matrix")
  mu_pop_draws <- fit$draws("mu_pop", format = "matrix")
  
  puma_summary_mrpr <- make_summary(mu_ps_draws,  PUMA_lev, "mrp-r")
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-p")
  
  return(list(
    puma_summary_mrpr=puma_summary_mrpr,
    puma_summary_mrpp=puma_summary_mrpp
  ))
  
}





getMRP_INT <- function(MR,
                      ps,
                      acs_pop,
                      mod = mod           
                      ) {
  
 
  grainsize = 500
  bin_fun = function(p) round(p, 1)
  psi_eps = 1e-6
  
   ## 0) LOCK LEVELS FROM TRAINING --------------------------------------------
  nps_ca <- MR
  nps_ca$AGEP  <- factor(nps_ca$AGEP)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev     <- levels(factor(nps_ca$PUMA))  
  
  # Fixed-effects design (no intercept; matches Stan)
  X_train <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train); n <- nrow(X_train)
  y <- as.integer(nps_ca$HICOV)
  
  
  group <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  J <- length(PUMA_lev)
  
  ## Helper to rebuild prediction matrices & group using training levels
  mk_pred <- function(df, X_train_cols, PUMA_lev, train_ref) {
    if (is.null(df) || nrow(df) == 0) {
      return(list(
        N = 0,
        X = matrix(0, 0, length(X_train_cols), dimnames = list(NULL, X_train_cols)),
        g = integer(0),
        df = df[0, , drop = FALSE]
      ))
    }
    df$AGEP  <- factor(df$AGEP,  levels = levels(train_ref$AGEP))
    df$SEX   <- factor(df$SEX,   levels = levels(train_ref$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(train_ref$RAC1P))
    df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
    
    Xp <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = df)
    
    miss <- setdiff(colnames(X_train_cols), colnames(Xp))
    if (length(miss)) {
      Xp <- cbind(Xp, matrix(0, nrow(Xp), length(miss),
                             dimnames = list(NULL, miss)))
    }
    Xp <- Xp[, X_train_cols, drop = FALSE]
    list(N = nrow(Xp), X = Xp, g = as.integer(df$PUMA), df = df)
  }
  
  ## 1) SELECTION MODEL (estimate ψ on stacked MR vs ps) ----------------------
  # Build a clean selection dataset with shared covariates
  sel_MR <- nps_ca[, c("AGEP", "SEX", "RAC1P", "PUMA")]
  sel_MR$S <- 1L
  sel_ps <- ps[, c("AGEP", "SEX", "RAC1P", "PUMA", "PWGTP")]
  sel_ps$S <- 0L
  
  # Align factor levels to training for selection model too
  align_levels <- function(df) {
    df$AGEP  <- factor(df$AGEP,  levels = levels(nps_ca$AGEP))
    df$SEX   <- factor(df$SEX,   levels = levels(nps_ca$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(nps_ca$RAC1P))
    df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
    df
  }
  sel_MR <- align_levels(sel_MR)
  sel_ps <- align_levels(sel_ps)
  
  sel_dat <- rbind(
    cbind(sel_MR, PWGTP = 1),
    sel_ps
  )
  
  sel_fit <- stats::glm(
    S ~ AGEP + SEX + RAC1P + PUMA,
    data = sel_dat,
    family = binomial(),
    weights = sel_dat$PWGTP
  )
  
  predict_psi <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(list(psi = numeric(0), psi_bin = integer(0)))
    df2 <- align_levels(df)
    p <- stats::predict(sel_fit, newdata = df2, type = "response")
    p <- pmin(pmax(p, psi_eps), 1 - psi_eps)
    list(psi = as.numeric(p), bin_val = bin_fun(p))
  }
  
  psi_train <- predict_psi(nps_ca)
  psi_ps    <- predict_psi(ps)
  psi_pop   <- predict_psi(acs_pop)
  
  # Build a consistent bin map across ALL datasets
  all_bin_vals <- unique(c(psi_train$bin_val, psi_ps$bin_val, psi_pop$bin_val))
  all_bin_vals <- sort(unique(all_bin_vals))
  bin_map <- setNames(seq_along(all_bin_vals), all_bin_vals)
  
  map_bins <- function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])
  
  psi_bin_train <- map_bins(psi_train$bin_val)
  psi_bin_ps    <- map_bins(psi_ps$bin_val)
  psi_bin_pop   <- map_bins(psi_pop$bin_val)
  G <- length(all_bin_vals)
  
  ## 2) BUILD Ps / Pop BLOCKS (design + groups + weights) ---------------------
  ps_b  <- mk_pred(ps,       colnames(X_train), PUMA_lev, nps_ca)
  p_b   <- mk_pred(acs_pop,  colnames(X_train), PUMA_lev, nps_ca)
  
  # Weights
  w_ps  <- if (ps_b$N > 0) as.vector(ps$PWGTP) else numeric(0)
  w_pop <- if (p_b$N > 0) rep(1, p_b$N) else numeric(0)
  
  ## 3) ASSEMBLE STAN DATA FOR MRP-INT ---------------------------------------
  stan_data <- list(
    # Training
    n = n,
    k = k,
    X = X_train,
    y = y,
    J = J,
    group = group,
    grainsize = grainsize,
    
    # psi for training (as probabilities; Stan takes logit() inside or we precompute there)
    psi = as.vector(psi_train$psi),
    G = G,
    psi_bin = psi_bin_train,
    
    # Ps (probability-sample poststratification; MRP-R summary)
    Nps = ps_b$N,
    Xps = ps_b$X,
    group_ps = ps_b$g,
    psi_ps = as.vector(psi_ps$psi),
    psi_bin_ps = psi_bin_ps,
    w_ps = w_ps,
    
    # Pop (full-population poststratification; MRP-P summary)
    Npop = p_b$N,
    Xpop = p_b$X,
    group_pop = p_b$g,
    psi_pop = as.vector(psi_pop$psi),
    psi_bin_pop = psi_bin_pop,
    w_pop = w_pop
  )
  
  ## 4) FIT STAN (MRP-INT) ----------------------------------------------------
  fit <- mod$sample(
    data = stan_data,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  )
  
  sum_tbl <- fit$summary()

  
  
  ## 5) COLLECT DRAWS & SUMMARIZE --------------------------------------------
  # mu_ps_draws  <- fit$draws("mu_ps",  format = "matrix")
  mu_pop_draws <- fit$draws("mu_pop", format = "matrix")
  
  make_summary <- function(draw_mat, group_names, tag) {
    draw_mat[!is.finite(draw_mat)] <- NA_real_
    stats_fn <- function(d) {
      qs <- stats::quantile(d, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      c(point_est = mean(d, na.rm = TRUE),
        sd        = stats::sd(d, na.rm = TRUE),
        lower_CI  = qs[1],
        upper_CI  = qs[2])
    }
    sm <- t(apply(draw_mat, 2, stats_fn))
    out <- tibble::tibble(PUMA = group_names) |>
      dplyr::bind_cols(as.data.frame(sm)) |>
      dplyr::mutate(model = tag, .before = 1)
    out
  }
  
  # puma_summary_mrpr <- make_summary(mu_ps_draws,  PUMA_lev, "mrp-int-ps")
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-int")
  
  list(
    fit = fit,
    # puma_summary_mrpr = puma_summary_mrpr,
    puma_summary_mrpp = puma_summary_mrpp,
    psi_bins = list(values = all_bin_vals, map = bin_map),
    rhat= sum_tbl$rhat
  )
}

