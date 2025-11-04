#' Multilevel Regression and Poststratification (MRP)
#'
#' Computes BOTH MRP-R and MRP-P estimates in a single Stan model run.
#'
#' @description
#' Fits a hierarchical logistic regression model to nonprobability sample data,
#' then post-stratifies predictions in TWO ways:
#'   1. MRP-R: Post-stratifies over probability sample (ps) cells, weighted by PWGTP
#'   2. MRP-P: Post-stratifies over full population (acs_pop) cells, equally weighted
#'
#' The Stan model (si2.stan) computes both estimates in its generated quantities block.
#'
#' @param MR Nonprobability sample data frame (e.g., nps) with outcome HICOV
#' @param ps Probability sample data frame for MRP-R post-stratification (must have PWGTP)
#' @param acs_pop Full population data frame for MRP-P post-stratification
#'
#' @return List with two data frames:
#'   - puma_summary_mrpr: MRP-R estimates (post-strat on ps)
#'   - puma_summary_mrpp: MRP-P estimates (post-strat on acs_pop)
#'
#' @details
#' Model: HICOV ~ AGEP + SEX + RAC1P + (1|PUMA)
#' Post-stratification cells: All unique combinations of (AGEP, SEX, RAC1P, PUMA)
#'                            from either ps (for MRP-R) or acs_pop (for MRP-P)
getMRP <- function(MR = nps,
                   ps = ps,
                   acs_pop = acs_pop) {
  ## 0) PREPARE TRAINING DATA --------------------------------------------
  nps_ca <- MR # nps

  # lock factor levels from training
  nps_ca$AGEP <- factor(nps_ca$AGEP)
  nps_ca$SEX <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev <- levels(factor(nps_ca$PUMA)) # lock PUMA levels

  # fixed-effects design
  X_train <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train)
  n <- nrow(X_train)
  y <- as.integer(nps_ca$HICOV)
  group <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  J <- length(PUMA_lev)

  ## 1) HELPER FUNCTION FOR PREDICTION MATRICES ------------------------------
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
    df$AGEP <- factor(df$AGEP, levels = levels(train_ref$AGEP))
    df$SEX <- factor(df$SEX, levels = levels(train_ref$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(train_ref$RAC1P))
    df$PUMA <- factor(df$PUMA, levels = PUMA_lev)

    Xp <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = df)
    # add any missing training columns (if some levels absent in this df)
    miss <- setdiff(X_train_cols, colnames(Xp))
    if (length(miss)) {
      Xp <- cbind(Xp, matrix(0, nrow(Xp), length(miss),
        dimnames = list(NULL, miss)
      ))
    }
    # order columns to match training
    Xp <- Xp[, X_train_cols, drop = FALSE]

    list(
      N = nrow(Xp),
      X = Xp,
      g = as.integer(df$PUMA),
      df_out = df
    )
  }

  ## 2) BUILD POST-STRATIFICATION CELLS --------------------------------------
  # This is where post-stratification cells are defined!
  # Cells = all unique combinations of (AGEP, SEX, RAC1P, PUMA) in each dataset

  # MRP-R: Post-stratification cells from probability sample
  # - Uses actual ps data rows as cells
  # - Each cell weighted by PWGTP (survey design weight)
  # - Accounts for sampling design but treats ps as fixed (no WFPBB)
  ps_b <- mk_pred(ps, colnames(X_train), PUMA_lev, nps_ca)
  w_ps <- ps$PWGTP[match(rownames(ps_b$X), rownames(ps))]

  # MRP-P: Post-stratification cells from full population
  # - Uses actual acs_pop data rows as cells
  # - Each cell weighted equally (w=1)
  # - Treats population as known/fixed
  p_b <- mk_pred(acs_pop, colnames(X_train), PUMA_lev, nps_ca)
  w_pop <- rep(1, p_b$N)

  ## 3) ASSEMBLE STAN DATA ---------------------------------------------------
  stan_data <- list(
    # training
    n = n,
    k = k,
    X = X_train,
    y = y,
    J = J,
    group = group,
    grainsize = 500 # choose 250–1000; tune for speed
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

  ## 4) FIT STAN MODEL -------------------------------------------------------
  fit <- mod$sample(
    data = stan_data,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123
  ) # 176.1

  ## 5) EXTRACT POST-STRATIFICATION RESULTS & SUMMARIZE ---------------------
  # NOTE: Post-stratification happens inside Stan's generated quantities block!
  # Stan computes weighted PUMA-level means using the cells we defined above:
  #   - mu_ps:  MRP-R estimates (aggregates predictions over ps cells with PWGTP weights)
  #   - mu_pop: MRP-P estimates (aggregates predictions over acs_pop cells equally)
  # These are already PUMA-level summaries, not individual predictions

  # Extract PUMA-level draws from Stan's generated quantities
  # mu_ps:  MRP-R PUMA means (post-stratified over ps with PWGTP weights)
  # mu_pop: MRP-P PUMA means (post-stratified over acs_pop with equal weights)
  p_ps_draws <- fit$draws("mu_ps", format = "matrix")   # S x Nps
  p_pop_draws <- fit$draws("mu_pop", format = "matrix") # S x Npop

  make_summary <- function(draw_mat, puma_names, tag) {
    draw_mat[!is.finite(draw_mat)] <- NA_real_
    tibble(PUMA = puma_names) |>
      bind_cols(as.data.frame(t(apply(draw_mat, 2, function(d) {
        qs <- quantile(d, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
        c(
          point_est = mean(d, na.rm = TRUE),
          sd = sd(d, na.rm = TRUE),
          lower_CI = qs[1],
          upper_CI = qs[2]
        )
      })))) |>
      mutate(model = tag, .before = 1)
  }

  puma_summary_mrpr <- make_summary(p_ps_draws, PUMA_lev, "mrp-r")
  puma_summary_mrpp <- make_summary(p_pop_draws, PUMA_lev, "mrp-p")

  return(list(
    puma_summary_mrpr = puma_summary_mrpr,
    puma_summary_mrpp = puma_summary_mrpp
  ))
}


#' MRP with Integrated Weighting (MRP-INT)
#'
#' MRP with selection bias correction via inclusion probability modeling.
#'
#' @description
#' Extends standard MRP by:
#'   1. Estimating inclusion probabilities (propensity to be in nonprob sample)
#'      using stacked nonprob + prob samples (following Valliant 2019)
#'   2. Including propensity scores in the outcome model as:
#'      - Linear term: beta_psi * logit(psi)
#'      - Categorical adjustment: random effects by propensity bin
#'   3. Post-stratifying over full population (acs_pop) cells
#'
#' This provides doubly-robust protection against model misspecification
#' (Si 2023, Section 2.3).
#'
#' @param MR Nonprobability sample data frame with outcome HICOV
#' @param ps Probability sample data frame (used to estimate inclusion probabilities)
#' @param acs_pop Full population data frame for post-stratification
#' @param mod Compiled Stan model object (mrp_int2.stan)
#'
#' @return List with:
#'   - puma_summary_mrpp: MRP-INT estimates (post-strat on acs_pop)
#'   - fit: Full Stan fit object
#'   - psi_bins: Propensity score binning information
#'   - rhat: Convergence diagnostics
#'
#' @details
#' Model: HICOV ~ AGEP + SEX + RAC1P + (1|PUMA) + beta_psi*logit(psi) + (1|psi_bin)
#' Inclusion probability model: P(in nonprob sample) ~ AGEP + SEX + RAC1P + PUMA
#' Post-stratification cells: All unique combinations in acs_pop
getMRP_INT <- function(MR,
                       ps,
                       acs_pop,
                       mod = mod) {
  grainsize <- 500
  bin_fun <- function(p) round(p, 1)
  psi_eps <- 1e-6

  ## 0) PREPARE TRAINING DATA --------------------------------------------
  # lock factor levels from training
  nps_ca <- MR
  nps_ca$AGEP <- factor(nps_ca$AGEP)
  nps_ca$SEX <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev <- levels(factor(nps_ca$PUMA))

  # Fixed-effects design (no intercept; matches Stan)
  X_train <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train)
  n <- nrow(X_train)
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
    df$AGEP <- factor(df$AGEP, levels = levels(train_ref$AGEP))
    df$SEX <- factor(df$SEX, levels = levels(train_ref$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(train_ref$RAC1P))
    df$PUMA <- factor(df$PUMA, levels = PUMA_lev)

    Xp <- model.matrix(~ 0 + AGEP + SEX + RAC1P, data = df)

    miss <- setdiff(colnames(X_train_cols), colnames(Xp))
    if (length(miss)) {
      Xp <- cbind(Xp, matrix(0, nrow(Xp), length(miss),
        dimnames = list(NULL, miss)
      ))
    }
    Xp <- Xp[, X_train_cols, drop = FALSE]
    list(N = nrow(Xp), X = Xp, g = as.integer(df$PUMA), df = df)
  }

  ## 1) ESTIMATE INCLUSION PROBABILITIES (Selection Model) -------------------
  # Following Valliant (2019) approach:
  # Step 1: Stack nonprob sample (S=1) with prob sample (S=0)
  # Step 2: Weight prob sample by PWGTP, nonprob sample by 1
  # Step 3: Fit weighted logistic regression P(S=1 | covariates)
  #
  # This estimates: P(unit is in nonprobability sample | AGEP, SEX, RAC1P, PUMA)
  # These are the "inclusion probabilities" or "propensity scores" (psi)

  sel_MR <- nps_ca[, c("AGEP", "SEX", "RAC1P", "PUMA")]
  sel_MR$S <- 1L # Nonprob sample indicator

  sel_ps <- ps[, c("AGEP", "SEX", "RAC1P", "PUMA", "PWGTP")]
  sel_ps$S <- 0L # Prob sample indicator

  # Align factor levels to training
  align_levels <- function(df) {
    df$AGEP <- factor(df$AGEP, levels = levels(nps_ca$AGEP))
    df$SEX <- factor(df$SEX, levels = levels(nps_ca$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(nps_ca$RAC1P))
    df$PUMA <- factor(df$PUMA, levels = PUMA_lev)
    df
  }
  sel_MR <- align_levels(sel_MR)
  sel_ps <- align_levels(sel_ps)

  # Stack: nonprob (weight=1) + prob (weight=PWGTP)
  sel_dat <- rbind(
    cbind(sel_MR, PWGTP = 1),
    sel_ps
  )

  # Fit inclusion probability model
  sel_fit <- stats::glm(
    S ~ AGEP + SEX + RAC1P + PUMA,
    data = sel_dat,
    family = binomial(),
    weights = sel_dat$PWGTP
  )

  predict_psi <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(list(psi = numeric(0), psi_bin = integer(0)))
    }
    df2 <- align_levels(df)
    p <- stats::predict(sel_fit, newdata = df2, type = "response")
    p <- pmin(pmax(p, psi_eps), 1 - psi_eps)
    list(psi = as.numeric(p), bin_val = bin_fun(p))
  }

  # Predict inclusion probabilities for all datasets
  psi_train <- predict_psi(nps_ca) # For training data
  psi_ps <- predict_psi(ps) # For probability sample (not used in final output)
  psi_pop <- predict_psi(acs_pop) # For population (used in MRP-INT post-strat)

  # Propensity score binning: Creates discrete categories for hierarchical modeling
  # - Bins propensities into ~11 groups: 0.0, 0.1, 0.2, ..., 1.0 (via round(p, 1))
  # - Allows Stan to estimate random effects per bin (zeta[psi_bin])
  # - More flexible than treating propensity as purely linear
  all_bin_vals <- unique(c(psi_train$bin_val, psi_ps$bin_val, psi_pop$bin_val))
  all_bin_vals <- sort(unique(all_bin_vals))
  bin_map <- setNames(seq_along(all_bin_vals), all_bin_vals)

  map_bins <- function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])

  psi_bin_train <- map_bins(psi_train$bin_val)
  psi_bin_ps <- map_bins(psi_ps$bin_val)
  psi_bin_pop <- map_bins(psi_pop$bin_val)
  G <- length(all_bin_vals) # Number of unique propensity bins

  ## 2) BUILD POST-STRATIFICATION CELLS (with propensity scores) -------------
  # Post-stratification cells = all unique combinations in acs_pop
  # Same as getMRP() but now each cell also has:
  #   - Estimated inclusion probability (psi)
  #   - Propensity bin assignment (psi_bin)
  # These are used as predictors in the outcome model

  ps_b <- mk_pred(ps, colnames(X_train), PUMA_lev, nps_ca)
  p_b <- mk_pred(acs_pop, colnames(X_train), PUMA_lev, nps_ca)

  # Weights (Note: MRP-INT only returns MRP-P style output, not MRP-R)
  w_ps <- if (ps_b$N > 0) as.vector(ps$PWGTP) else numeric(0)
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
      c(
        point_est = mean(d, na.rm = TRUE),
        sd = stats::sd(d, na.rm = TRUE),
        lower_CI = qs[1],
        upper_CI = qs[2]
      )
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
    rhat = sum_tbl$rhat
  )
}
