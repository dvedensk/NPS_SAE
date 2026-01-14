
## MRP
# ---------- 1) individuals to cells ----------
.align_to_train <- function(df, train_ref, PUMA_lev) {
  df$AGEP_binned  <- factor(df$AGEP_binned,  levels = levels(train_ref$AGEP_binned))
  df$SEX   <- factor(df$SEX,   levels = levels(train_ref$SEX))
  df$RAC1P <- factor(df$RAC1P, levels = levels(train_ref$RAC1P))
  df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
  df
}

collapse_to_cells <- function(df, train_ref, PUMA_lev,
                              design_formula = ~ 0 + AGEP_binned + SEX + RAC1P,
                              weight_col = NULL) {
  df <- .align_to_train(df, train_ref, PUMA_lev)
  w <- if (!is.null(weight_col) && weight_col %in% names(df)) df[[weight_col]] else 1
  w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
  
  cells <- df |>
    dplyr::mutate(W = w) |>
    dplyr::group_by(AGEP_binned, SEX, RAC1P, PUMA) |>
    dplyr::summarise(w = sum(W), .groups = "drop")
  
  Xp <- model.matrix(design_formula, data = cells)
  g  <- as.integer(cells$PUMA)                # 1..J
  list(Xp = Xp, g = g, w = cells$w, cells = cells)
}

# ---------- 2) extract CmdStanR draws ----------
get_beta_draws <- function(fit) {
  # S x K (format='matrix' gives draws as rows)
  fit$draws("beta", format = "matrix")
}

get_apuma_draws <- function(fit) {
  ap <- fit$draws("a_puma", format = "matrix")  # S x J (as columns a_puma[1]...)
  ap
}
# ---------- 3) poststrat by cell ----------
poststrat_SxJ <- function(beta, a_puma,Xp, g, w = NULL, J = max(g)) {
  # beta: S x K   (posterior draws)
  # Xp  : N x K   (prediction design matrix)
  # g   : length-N integer group indices 1..J
  # w   : length-N weights (if NULL, all 1)
  # J   : number of groups (PUMAs)
  
  
  N <- nrow(Xp)
  S <- nrow(beta)
  
  # weights
  if (is.null(w)) w <- rep(1, N) else w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  # 1) Compute probabilities for *all cells × all draws*
  #    result: (N x S)
  eta <- Xp %*% t(beta)
  # Add random intercept per cell: for each cell i, add a_puma[, g[i]] (length S)
  # a_puma[, g] gives S x N, transpose -> N x S
  eta <- eta + t(a_puma[, g, drop = FALSE])
  
  p   <- plogis(eta)# N x S
  
  # 2) Pre-split indices and normalize weights by PUMA
  idx_by_j <- split(seq_len(N), factor(g, levels = seq_len(J)))
  w_norm <- lapply(idx_by_j, function(idx) {
    if (length(idx) == 0) return(numeric(0))
    sw <- sum(w[idx])
    if (sw > 0) w[idx] / sw else rep(0, length(idx))
  })
  
  # 3) Aggregate weighted means for each PUMA (S x J)
  mu <- matrix(NA_real_, nrow = S, ncol = J)
  for (j in seq_len(J)) {
    idx <- idx_by_j[[j]]
    if (length(idx) == 0) next
    wj <- w_norm[[j]]
    mu[, j] <- as.numeric(crossprod(wj, p[idx, , drop = FALSE]))
  }
  
  mu  # (S x J)
}

# ---------- 4) summary ----------

make_summary <- function(draw_mat, puma_names, tag) {
  draw_mat[!is.finite(draw_mat)] <- NA_real_
  tibble::tibble(PUMA = puma_names) |>
    dplyr::bind_cols(as.data.frame(t(apply(draw_mat, 2, function(d) {
      qs <- stats::quantile(d, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      c(point_est = mean(d, na.rm = TRUE),
        sd        = stats::sd(d,   na.rm = TRUE),
        lower_CI  = qs[1],
        upper_CI  = qs[2])
    })))) |>
    dplyr::mutate(model = tag, .before = 1)
}


# ---------- 5) getMRP---------- 

getMRP=function(MR=nps,
                ps=ps,
                acs_pop=acs_pop,
                bootstrap=FALSE,
                L=5,
                threads=4){

  # TRAINING DATA ----
  nps_ca <- MR  # nps
  
  # lock factor levels from training
  nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev <- levels(factor(nps_ca$PUMA))  # lock PUMA levels
  
  # fixed-effects design 
  X_train <- model.matrix(~ 0 + AGEP_binned + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train)
  n <- nrow(X_train)
  y <- as.integer(nps_ca$PUBCOV)


  #threads <- 4
  grainsize=as.integer(n / (threads))

  #  data dependent prior for beta
  sd_y <- stats::sd(as.numeric(y))
  sd_x <- apply(X_train, 2, stats::sd)
  sd_x[sd_x == 0] <- NA_real_   # or small number
  beta_scale <- 2.5 * sd_y / sd_x
  
  
  puma_id <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  J <- length(PUMA_lev)
  
  stan_data <- list(
    n = n,
    k = k,
    X = X_train,
    y = y,
    J = J,
    puma_id = puma_id,
    grainsize = grainsize,
    beta_scale = as.vector(beta_scale)
  )
  

  
  fit <- mod$sample(
    data = stan_data,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    seed = 123
  )
  
  sum_tbl <- fit$summary()
  
  
  # fit_glm <- rstanarm::stan_glm( #stan_glmer
  #   PUBCOV ~ 0 + AGEP_binned + SEX + RAC1P+ WAGP+ACCESSINET, #(1|PUMA)
  #   data   = nps_ca,
  #   family = binomial(link = "logit"),
  #   chains = 2,
  #   iter   = 1500,
  #   warmup = 500,
  #   cores  = 2,
  #   seed   = 123
  # )

  
  # TRAINING ref (to lock levels & PUMA order)
  nps_ca <- MR
  nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev     <- levels(factor(nps_ca$PUMA))
  J            <- length(PUMA_lev)
  
  # Collapse Ps and Pop to cells (huge speed/memory win)
  ps_cells  <- collapse_to_cells(ps,      nps_ca, PUMA_lev, weight_col = "weights")
  pop_cells <- collapse_to_cells(acs_pop, nps_ca, PUMA_lev, weight_col = "weights")
  
  beta <- get_beta_draws(fit)   
  apuma_draws <- get_apuma_draws(fit)
  
  mu_ps_draws  <- if (nrow(ps_cells$Xp)  > 0)
    poststrat_SxJ(beta, apuma_draws, ps_cells$Xp,  ps_cells$g,  ps_cells$w,  J) else
      matrix(NA_real_, nrow = nrow(beta), ncol = J)
  
  mu_pop_draws <- if (nrow(pop_cells$Xp) > 0)
    poststrat_SxJ(beta=beta, 
                  apuma_draws,
                  Xp=pop_cells$Xp, 
                  g=pop_cells$g, 
                  w=pop_cells$w, 
                  J=J) else
                    matrix(NA_real_, nrow = nrow(beta), ncol = J)
  

  puma_summary_mrpr <- make_summary(mu_ps_draws,  PUMA_lev, "mrp-r")
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-p")


  if(bootstrap){
    # Preallocate array for S × J × L (posterior × PUMAs × bootstrap)
    S <- nrow(beta)
    mu_boot_array <- array(NA_real_, dim = c(S, J, L))

    cat("\nGenerating", L, "bootstrap samples...\n")
    system.time({
      for (l in 1:L) {
        # Generate bootstrap sample l using simple weighted bootstrap
        # (much faster than MSIMST::WFPBB which fits a Bayesian model)
        boot_idx <- sample(
          1:nrow(ps),
          size = nrow(ps),
          replace = TRUE,
          prob = ps$weights
        )

        # Resample THIS bootstrap replicate
        ps_boot <- ps[boot_idx, ]

        # Collapse THIS bootstrap sample to cells
        ps_boot_cells <- collapse_to_cells(ps_boot, nps_ca, PUMA_lev, weight_col = "weights")

        # Poststratify THIS bootstrap sample
        mu_boot_array[, , l] <- if (nrow(ps_boot_cells$Xp) > 0)
          poststrat_SxJ(beta=beta,
                        apuma_draws,
                        Xp=ps_boot_cells$Xp,
                        g=ps_boot_cells$g,
                        w=ps_boot_cells$w,
                        J=J) else
                          matrix(NA_real_, nrow = S, ncol = J)

        if (l %% 10 == 0) cat("  Completed", l, "of", L, "bootstrap samples\n")
      }
    })

    # Flatten S × J × L into (S*L) × J to capture both uncertainties
    mu_combined <- matrix(aperm(mu_boot_array, c(3, 1, 2)), nrow = S*L, ncol = J)
    puma_summary_mrpr_bootstrap <- make_summary(mu_combined,  PUMA_lev, "mrp-r-bootstrap")


  }else{
    puma_summary_mrpr_bootstrap=NULL
  }
  
  
  return(list(
    puma_summary_mrpr=puma_summary_mrpr,
    puma_summary_mrpp=puma_summary_mrpp,
    puma_summary_mrpr_bootstrap=puma_summary_mrpr_bootstrap,
    rhat=sum_tbl$rhat
  ))
  
}

## MRP_INT
# ---------- 0) tiny utilites ----------  


bin_fun <- function(p, digits = 2) round(p, digits = digits)
# map_bins <- function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])

#round(p, 1)  may end up only has two groups if inclusion probs are small
.align_to_train <- function(df, train_ref, PUMA_lev) {
  df$AGEP_binned  <- factor(df$AGEP_binned,  levels = levels(train_ref$AGEP_binned))
  df$SEX   <- factor(df$SEX,   levels = levels(train_ref$SEX))
  df$RAC1P <- factor(df$RAC1P, levels = levels(train_ref$RAC1P))
  df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
  df
}

# ---------- 1) individuals to cells ----------

collapse_to_cells_int <- function(df, train_ref, PUMA_lev,
                                  psi_vec,
                                  psi_bin_vec,
                                  design_formula = ~ 0 + AGEP_binned + SEX + RAC1P,
                                  weight_col = NULL,
                                  bin_map = bin_map,
                                  map_bins = function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])) {
  
  df <- .align_to_train(df, train_ref, PUMA_lev)
  
  w <- if (!is.null(weight_col) && weight_col %in% names(df)) df[[weight_col]] else 1
  w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
  
  if (length(psi_vec) != nrow(df) || length(psi_bin_vec) != nrow(df)) {
    stop("psi_vec and psi_bin_vec must have length nrow(df)")
  }
  df$psi     <- as.numeric(psi_vec)
  df$psi_bin <- as.integer(psi_bin_vec)
  
  # critical: PUMA levels locked
  df$PUMA <- factor(df$PUMA, levels = PUMA_lev)
  
  cells <- df |>
    dplyr::mutate(W = w) |>
    dplyr::group_by(AGEP_binned, SEX, RAC1P, PUMA) |>
    dplyr::summarise(
      w   = sum(W),
      psi = mean(psi),
      .groups = "drop"
    )
  
  cells$PUMA <- factor(cells$PUMA, levels = PUMA_lev)
  
  # recompute bins from cell-mean psi (your logic)
  cells$psi_bin <- map_bins(bin_fun(cells$psi,digits = 2))
  
  Xp <- model.matrix(design_formula, data = cells)
  g  <- as.integer(cells$PUMA)
  if (anyNA(g)) stop("Some cells have PUMA not in PUMA_lev -> NA group id.")
  
  list(Xp = Xp, g = g, w = cells$w, psi = cells$psi, psi_bin = cells$psi_bin, cells = cells)
}

# ---------- 2) Extract CmdStanR draws ----------
get_params_int_draws <- function(fit) {
  list(
    beta     = fit$draws("beta", format = "matrix"),      # S x K
    beta_psi = as.numeric(fit$draws("beta_psi", format = "matrix")),  # S
    zeta     = fit$draws("zeta", format = "matrix"),      # S x G
    a_puma   = fit$draws("a_puma", format = "matrix")     # S x J
  )
}

# ---------- 3) poststrat by cell ----------
poststrat_int_SxJ <- function(beta, beta_psi, zeta, a_puma,
                              Xp, g, psi, psi_bin, w = NULL, J = max(g)) {
  N <- nrow(Xp)
  S <- nrow(beta)
  
  if (is.null(w)) w <- rep(1, N) else w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  lp <- qlogis(pmin(pmax(psi, 1e-6), 1 - 1e-6))  # N
  
  # N x S components
  eta <- Xp %*% t(beta)                                  # fixed
  eta <- eta + tcrossprod(lp, beta_psi)                  # selection slope
  eta <- eta + t(zeta[, psi_bin, drop = FALSE])          # psi-bin 
  eta <- eta + t(a_puma[, g, drop = FALSE])              # PUMA
  
  p <- plogis(eta)  # N x S
  
  idx_by_j <- split(seq_len(N), factor(g, levels = seq_len(J)))
  mu <- matrix(NA_real_, nrow = S, ncol = J)
  
  for (j in seq_len(J)) {
    idx <- idx_by_j[[j]]
    if (length(idx) == 0) next
    sw <- sum(w[idx])
    wj <- if (sw > 0) w[idx] / sw else rep(0, length(idx))
    mu[, j] <- as.numeric(crossprod(wj, p[idx, , drop = FALSE]))
  }
  mu
}



# ---------- 4) getMRP_INT---------- 
getMRP_INT <- function(MR,
                       ps,
                       acs_pop,
                       mod = mod,
                       adjust=TRUE,
                       bootstrap=FALSE,
                       L=5,
                       threads=4
) {
  
  
  # ps$PWGTP<-ps$weights
  
  
  psi_eps = 1e-6
  
  ## 0) LOCK LEVELS FROM TRAINING --------------------------------------------
  nps_ca <- MR
  nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev     <- levels(factor(nps_ca$PUMA))  
  J            <- length(PUMA_lev) 
  
  # Fixed-effects design (no intercept; matches Stan)
  X_train <- model.matrix(~ 0 + AGEP_binned + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train)
  n <- nrow(X_train)
  y <- as.integer(nps_ca$PUBCOV)
  
  
  
  ## 1) SELECTION MODEL --------------------------------------------
  # Build a clean selection dataset with shared covariates
  sel_MR <- nps_ca[, c("AGEP_binned", "SEX", "RAC1P", "PUMA")]
  sel_MR$S <- 1L
  sel_ps <- ps[, c("AGEP_binned", "SEX", "RAC1P", "PUMA", "weights")]
  sel_ps$S <- 0L
  
  # Align factor levels to training for selection model too
  align_levels <- function(df) {
    df$AGEP  <- factor(df$AGEP_binned,  levels = levels(nps_ca$AGEP_binned))
    df$SEX   <- factor(df$SEX,   levels = levels(nps_ca$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(nps_ca$RAC1P))
    df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
    df
  }
  sel_MR <- align_levels(sel_MR)
  sel_ps <- align_levels(sel_ps)
  
  sel_dat <- rbind(
    cbind(sel_MR, weights = 1),
    sel_ps
  )
  
  if(adjust==TRUE){
    ## 1) Estimated population size from PS
    N_hat <- sum(sel_ps$weights)
    ## 2) Size of nonprobability sample
    n_np <- nrow(sel_MR)
    ## 3) Adjustment factor
    adj_factor <- (N_hat - n_np) / N_hat
    ## 4) Adjust PS weights
    sel_ps$weights <- sel_ps$weights * adj_factor
    sel_MR$weights=1
    sel_dat <- rbind(
      sel_MR[, c("AGEP_binned","SEX","RAC1P","PUMA","S","weights")],
      sel_ps[, c("AGEP_binned","SEX","RAC1P","PUMA","S","weights")]
    ) 
  }else{
    adj_factor=NULL
  }
  
  sel_fit <- stats::glm(
    S ~ AGEP_binned + SEX + RAC1P, 
    data = sel_dat,
    family = binomial(),
    weights = sel_dat$weights
  )
  
  #bin_val is 0.11, 0.12, etc, psi is the predicted psi, using bin_val represents the random-effect group which is used during poststratification.
  predict_psi <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(list(psi = numeric(0), psi_bin = integer(0)))
    df2 <- align_levels(df)
    p <- stats::predict(sel_fit, newdata = df2, type = "response")
    p <- pmin(pmax(p, psi_eps), 1 - psi_eps)#to avoid 0
    list(psi = as.numeric(p), bin_val = bin_fun(p,digits = 2))
  }
  
  psi_train <- predict_psi(nps_ca)
  psi_ps    <- predict_psi(ps)
  psi_pop   <- predict_psi(acs_pop)
  
  # Build a consistent bin map across ALL datasets
  all_bin_vals <- unique(c(psi_train$bin_val, psi_ps$bin_val, psi_pop$bin_val))
  all_bin_vals <- sort(unique(all_bin_vals))
  G <- length(all_bin_vals)
  
  
  bin_map <- setNames(seq_along(all_bin_vals), all_bin_vals)
  map_bins <- function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])

  
  psi_bin_train <- map_bins(psi_train$bin_val)
  psi_bin_ps    <- map_bins(psi_ps$bin_val)
  psi_bin_pop   <- map_bins(psi_pop$bin_val)

  
  
  
  pop_cells <- collapse_to_cells_int(
    df          = acs_pop,
    train_ref   = nps_ca,
    PUMA_lev    = PUMA_lev,
    # psi_bin_val = psi_pop$bin_val, # 0.3
    psi_vec     = psi_pop$psi,
    psi_bin_vec = psi_bin_pop, # 3
    weight_col  = "weights",
    bin_map = bin_map,
    map_bins= function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])
    
  )
  
  ps_cells <- collapse_to_cells_int(
    df          = ps,
    train_ref   = nps_ca,
    PUMA_lev    = PUMA_lev,
    psi_vec     = psi_ps$psi,
    psi_bin_vec = psi_bin_ps, # 3
    weight_col  = "weights",
    bin_map = bin_map,
    map_bins= function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])
    
  )
  
  
  ## 2) MAKE STAN DATA FOR MRP-INT ---------------------------------------
  # data-dependent prior scales
  
  # outcome sd (binary y)
  sd_y <- stats::sd(as.numeric(y)) 
  
  # column-wise sd for X
  sd_x <- apply(X_train, 2, stats::sd)
  sd_x[sd_x == 0] <- NA_real_   
  beta_scale <- 1 * (sd_y)  / sd_x 
  # when 2.5 * sd_y / sd_x is used, warning:  Exception: bernoulli_logit_glm_lpmf: Intercept[1] is -inf, but must be finite!
  # beta_scale[!is.finite(beta_scale)] <- 5
  # beta_scale <- pmin(beta_scale, 5) # to avoid intercept[1]=inf
  # 
  beta_psi_scale <- 2.5 * (sd_y) / stats::sd(qlogis(psi_train$psi))
  sigma_psi_rate <- 1 / sd_y   
  
  
  n <- nrow(X_train)
  # threads <- 4
  grainsize=as.integer(n / (threads*4))
  
  puma_id_train <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  J <- length(PUMA_lev)
  
  
  stan_data <- list(
    n = n,
    k = k,
    X = X_train,
    y = y,
    grainsize = grainsize,
    lp_psi = qlogis(psi_train$psi),                  
    G = G,
    psi_bin = psi_bin_train,
    
    
    J = J,
    puma_id = puma_id_train,
    
    beta_scale = as.vector(beta_scale),
    beta_psi_scale = beta_psi_scale,
    sigma_psi_rate = sigma_psi_rate
    
    # sigma_puma_rate = 1 / sd_y,
    
    
    
  )
  ## 3) FIT STAN (MRP-INT) ----------------------------------------------------
  fit <- mod$sample(
    data = stan_data,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = threads,
    iter_warmup = 500,
    iter_sampling = 1000
    # seed = 1233
  )
  
  
  sum_tbl <- fit$summary()

  
  ## 4) COLLECT DRAWS & SUMMARIZE --------------------------------------------
  
  # nps_ca <- MR
  # nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  # nps_ca$SEX   <- factor(nps_ca$SEX)
  # nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  # PUMA_lev     <- levels(factor(nps_ca$PUMA))
  # J            <- length(PUMA_lev)
  params <- get_params_int_draws(fit)
  beta     <- params$beta
  beta_psi <- params$beta_psi
  zeta     <- params$zeta
  a_puma   <- params$a_puma
  
  
  
  
  # Pop poststrat (MRP-P with INT correction)
  mu_pop_draws <- if (nrow(pop_cells$Xp) > 0) {
    poststrat_int_SxJ(
      beta = beta, 
      beta_psi = beta_psi, 
      zeta = zeta, 
      a_puma = a_puma,
      Xp = pop_cells$Xp, 
      g = pop_cells$g, 
      psi = pop_cells$psi,
      psi_bin = pop_cells$psi_bin,
      w = pop_cells$w,
      J = J
    )
  } else {
    matrix(NA_real_, nrow = nrow(beta), ncol = J)
  }
  
  
  mu_ps_draws <- if (nrow(pop_cells$Xp) > 0) {
    poststrat_int_SxJ(
      beta = beta, 
      beta_psi = beta_psi, 
      zeta = zeta, 
      a_puma = a_puma,
      Xp = ps_cells$Xp, 
      g = ps_cells$g, 
      psi = ps_cells$psi,
      psi_bin = ps_cells$psi_bin,
      w = ps_cells$w,
      J = J
    )
  } else {
    matrix(NA_real_, nrow = nrow(beta), ncol = J)
  }
  
  
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-p-int")
  puma_summary_mrpr <- make_summary(mu_ps_draws, PUMA_lev, "mrp-r-int")

  if(bootstrap){
    # Preallocate array for S × J × L (posterior × PUMAs × bootstrap)
    S <- nrow(beta)
    mu_boot_array <- array(NA_real_, dim = c(S, J, L))

    cat("\nGenerating", L, "bootstrap samples for MRP-INT...\n")
    system.time({
      for (l in 1:L) {
        # Generate bootstrap sample l using simple weighted bootstrap
        boot_idx <- sample(
          1:nrow(ps),
          size = nrow(ps),
          replace = TRUE,
          prob = ps$weights
        )

        # Resample THIS bootstrap replicate
        ps_boot <- ps[boot_idx, ]

        # Compute psi for bootstrap sample
        psi_ps_boot <- predict_psi(ps_boot)
        psi_bin_ps_boot <- map_bins(psi_ps_boot$bin_val)

        # Collapse THIS bootstrap sample to cells
        ps_boot_cells <- collapse_to_cells_int(
          df          = ps_boot,
          train_ref   = nps_ca,
          PUMA_lev    = PUMA_lev,
          psi_vec     = psi_ps_boot$psi,
          psi_bin_vec = psi_bin_ps_boot,
          weight_col  = "weights",
          bin_map     = bin_map,
          map_bins    = function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])
        )

        # Poststratify THIS bootstrap sample
        mu_boot_array[, , l] <- if (nrow(ps_boot_cells$Xp) > 0)
          poststrat_int_SxJ(
            beta = beta,
            beta_psi = beta_psi,
            zeta = zeta,
            a_puma = a_puma,
            Xp = ps_boot_cells$Xp,
            g = ps_boot_cells$g,
            psi = ps_boot_cells$psi,
            psi_bin = ps_boot_cells$psi_bin,
            w = ps_boot_cells$w,
            J = J
          ) else
            matrix(NA_real_, nrow = S, ncol = J)

        if (l %% 10 == 0) cat("  Completed", l, "of", L, "bootstrap samples\n")
      }
    })

    # Flatten S × J × L into (S*L) × J to capture both uncertainties
    mu_combined <- matrix(aperm(mu_boot_array, c(3, 1, 2)), nrow = S*L, ncol = J)
    puma_summary_mrpr_bootstrap <- make_summary(mu_combined, PUMA_lev, "mrp-r-int-bootstrap")

  }else{
    puma_summary_mrpr_bootstrap=NULL
  }

  list(
    fit          = fit,
    puma_summary_mrpp = puma_summary_mrpp,
    puma_summary_mrpr = puma_summary_mrpr,
    puma_summary_mrpr_bootstrap = puma_summary_mrpr_bootstrap,
    psi_train     = psi_train,
    rhat         = sum_tbl$rhat,
    adj_factor   = adj_factor
  )
}

