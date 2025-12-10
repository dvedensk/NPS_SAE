
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

# ---------- 3) poststrat by cell ----------
poststrat_SxJ <- function(beta, Xp, g, w = NULL, J = max(g)) {
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
  p   <- plogis(eta)
  
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
# Load compiled Stan models
mod <- cmdstan_model(
  file.path("/Users/qianyudong/Dropbox/simulation/code", "si2_1111.stan"),
  cpp_options = list(stan_threads = TRUE))



getMRP=function(MR=nps,
                ps=ps,
                acs_pop=acs_pop,
                WFPBB=FALSE,
                L=5){
  
  # TRAINING DATA ----
  nps_ca <- MR  # nps
  
  # lock factor levels from training
  nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev <- levels(factor(nps_ca$PUMA))  # lock PUMA levels
  
  # fixed-effects design 
  X_train <- model.matrix(~ 0 + AGEP_binned + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train); n <- nrow(X_train)
  y <- as.integer(nps_ca$PUBCOV)

  

  n <- nrow(X_train)
  threads <- 4
  grainsize=as.integer(n / (threads))

  #  data dependent prior for beta
  sd_y <- stats::sd(as.numeric(y))
  sd_x <- apply(X_train, 2, stats::sd)
  sd_x[sd_x == 0] <- NA_real_   # or small number
  beta_scale <- 2.5 * sd_y / sd_x
  
  stan_data <- list(
    # training
    n = n,
    k = k,
    X = X_train,
    y = y,
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
  # fit_glm <- rstanarm::stan_glm(
  #   HICOV ~ 0 + AGEP + SEX + RAC1P,    
  #   data   = nps_ca,
  #   family = binomial(link = "logit"),
  #   chains = 2,
  #   iter   = 1500,  
  #   warmup = 500,
  #   cores  = 2,
  #   seed   = 123
  # )
  # takes 72s. 
  

  

  
  
  # 
  # # Build Ps (probability-sample poststrat) block
  # ps_b <- mk_pred(ps, colnames(X_train), PUMA_lev, nps_ca)  
  # w_ps <- ps$PWGTP[match(rownames(ps_b$X), rownames(ps))]  
  # 
  # # Build P (full population poststrat) block
  # p_b  <- mk_pred(acs_pop,  colnames(X_train), PUMA_lev, nps_ca)
  # w_pop <- rep(1, p_b$N)
  # 
  
  # TRAINING ref (to lock levels & PUMA order)
  nps_ca <- MR
  nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev     <- levels(factor(nps_ca$PUMA))
  J            <- length(PUMA_lev)
  
  # Collapse Ps and Pop to cells (huge speed/memory win)
  ps_cells  <- collapse_to_cells(ps,      nps_ca, PUMA_lev, weight_col = "PWGTP")
  pop_cells <- collapse_to_cells(acs_pop, nps_ca, PUMA_lev, weight_col = "PWGTP")
  
  beta <- get_beta_draws(fit)   # S x K
  
  mu_ps_draws  <- if (nrow(ps_cells$Xp)  > 0)
    poststrat_SxJ(beta, ps_cells$Xp,  ps_cells$g,  ps_cells$w,  J) else
      matrix(NA_real_, nrow = nrow(beta), ncol = J)
  
  mu_pop_draws <- if (nrow(pop_cells$Xp) > 0)
    poststrat_SxJ(beta=beta, 
                  Xp=pop_cells$Xp, 
                  g=pop_cells$g, 
                  w=pop_cells$w, 
                  J=J) else
      matrix(NA_real_, nrow = nrow(beta), ncol = J)
  

  puma_summary_mrpr <- make_summary(mu_ps_draws,  PUMA_lev, "mrp-r")
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-p")
  
  
  if(WFPBB){
    index_list <- vector("list", L)
    system.time({
      for (l in 1:L) {
        index_list[[l]] <- MSIMST::WFPBB(
          y = 1:nrow(ps),
          w = ps$PWGTP,
          N = sum(ps$PWGTP),
          n = nrow(ps),
          verbatim = FALSE
        )
        print(paste("Now generate WFPBB", l))
      }
    })
    ps_star <- ps[unlist(index_list), ]
    ps_star_cells <- collapse_to_cells(ps_star, nps_ca, PUMA_lev, weight_col = "PWGTP")
    
    # Build Ps (probability-sample poststrat) block
    ps_b_star <- mk_pred(ps_star, colnames(X_train), PUMA_lev, nps_ca)  
    w_ps_star <- ps_star$PWGTP[match(rownames(ps_b_star$X), rownames(ps_star))]  
    
    mu_star_draws <- if (nrow(ps_star_cells$Xp) > 0)
      poststrat_SxJ(beta=beta, 
                    Xp=ps_star_cells$Xp, 
                    g=ps_star_cells$g, 
                    w=ps_star_cells$w, 
                    J=J) else
                      matrix(NA_real_, nrow = nrow(beta), ncol = J)
    
    puma_summary_mrpr_WFPBB <- make_summary(mu_star_draws,  PUMA_lev, "mrp-r-WFPBB")
    
    
  }else{
    puma_summary_mrpr_WFPBB=NULL
  }
  
  
  return(list(
    puma_summary_mrpr=puma_summary_mrpr,
    puma_summary_mrpp=puma_summary_mrpp,
    puma_summary_mrpr_WFPBB=puma_summary_mrpr_WFPBB,
    rhat=sum_tbl$rhat
  ))
  
}

## MRP_INT
# ---------- 0) tiny utilites ----------  


bin_fun <- function(p) round(p, 1)  
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
                                  bin_map= bin_map,
                                  map_bins= function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])) {
  # Align to training factors
  df <- .align_to_train(df, train_ref, PUMA_lev)
  
  # Weights
  w <- if (!is.null(weight_col) && weight_col %in% names(df)) df[[weight_col]] else 1
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  # Attach psi and psi_bin to rows
  if (length(psi_vec) != nrow(df) || length(psi_bin_vec) != nrow(df)) {
    stop("psi_vec and psi_bin_vec must have length nrow(df)")
  }
  df$psi     <- as.numeric(psi_vec)
  df$psi_bin <- as.integer(psi_bin_vec)
  
  
  cells <- df |>
    dplyr::mutate(W = w) |>
    dplyr::group_by(AGEP_binned, SEX, RAC1P, PUMA) |>
    dplyr::summarise(
      w       = sum(W),
      psi     = mean(psi),
      .groups = "drop"
    )
  cells$psi_bin<-bin_fun(cells$psi)
  cells$psi_bin<-map_bins( cells$psi_bin)
  
  
  
  
  # Design matrix for outcome model
  Xp <- model.matrix(design_formula, data = cells)
  g  <- as.integer(cells$PUMA)  # group index 1..J
  
  
  list(
    Xp    = Xp,
    g     = g,
    w     = cells$w,
    psi   = cells$psi,
    psi_bin = cells$psi_bin,
    cells = cells
  )
}
# ---------- 2) Extract CmdStanR draws ----------
get_params_int_draws <- function(fit) {
  # beta: S x K
  beta_mat <- fit$draws("beta", format = "matrix")
  
  # beta_psi: S-length vector
  beta_psi_mat <- fit$draws("beta_psi", format = "matrix")
  beta_psi_vec <- as.numeric(beta_psi_mat[, 1])
  
  # zeta: S x G   (columns like zeta[1], zeta[2], ...)
  zeta_mat <- fit$draws("zeta", format = "matrix")
  
  list(
    beta     = beta_mat,
    beta_psi = beta_psi_vec,
    zeta     = zeta_mat
  )
}

# ---------- 3) poststrat by cell ----------

poststrat_int_SxJ <- function(beta,
                              beta_psi,
                              zeta,
                              Xp,
                              g,
                              psi,
                              psi_bin,
                              w = NULL,
                              J = max(g)) {
  # beta    : S x K
  # beta_psi: length-S
  # zeta    : S x G
  # Xp      : N x K
  # g       : length-N group indices 1..J (PUMAs)
  # psi     : length-N probabilities in (0,1)
  # psi_bin : length-N integers in 1..G
  # w       : length-N weights
  # returns : S x J matrix of PUMA-level poststratified means
  
  N <- nrow(Xp)
  S <- nrow(beta)
  
  if (length(psi) != N) stop("psi must have length N = nrow(Xp)")
  if (length(psi_bin) != N) stop("psi_bin must have length N = nrow(Xp)")
  if (length(beta_psi) != S) stop("beta_psi must have length S = nrow(beta)")
  
  # weights
  if (is.null(w)) w <- rep(1, N) else w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  
  # logit(psi) for cells
  lp_psi <- qlogis(psi)
  
  # Fixed-effects part: N x S
  eta_fixed <- Xp %*% t(beta)   # (N x K) %*% (K x S) = N x S
  
  # Selection term: beta_psi_s * lp_psi_i
  eta_sel <- outer(lp_psi, beta_psi, "*")  # N x S
  
  # Psi-bin random effects: zeta_s,psi_bin[i]
  # zeta is S x G; transpose to G x S and index by psi_bin
  zeta_t   <- t(zeta)                              # G x S
  eta_psiRE <- zeta_t[psi_bin, , drop = FALSE]     # N x S
  
  # Total linear predictor and probabilities
  eta <- eta_fixed + eta_sel + eta_psiRE
  p   <- plogis(eta)  # N x S
  
  # Split indices by PUMA
  idx_by_j <- split(seq_len(N), factor(g, levels = seq_len(J)))
  
  # Normalized weights within PUMA
  w_norm <- lapply(idx_by_j, function(idx) {
    if (length(idx) == 0) return(numeric(0))
    sw <- sum(w[idx])
    if (sw > 0) w[idx] / sw else rep(0, length(idx))
  })
  
  # Aggregate S x J
  mu <- matrix(NA_real_, nrow = S, ncol = J)
  for (j in seq_len(J)) {
    idx <- idx_by_j[[j]]
    if (length(idx) == 0) next
    wj <- w_norm[[j]]
    # p[idx, ] is (length(idx) x S); crossprod gives S x 1
    mu[, j] <- as.numeric(crossprod(wj, p[idx, , drop = FALSE]))
  }
  
  mu
}

mod <- cmdstan_model(
  file.path("/Users/qianyudong/Dropbox/simulation/code", "mrp_int2_1111.stan"),
  cpp_options = list(stan_threads = TRUE),
  force_recompile = TRUE
)

# ---------- 4) getMRP_INT---------- 
getMRP_INT <- function(MR,
                       ps,
                       acs_pop,
                       mod = mod,
                       adjust=T
) {
  
  
  
  
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
  sel_ps <- ps[, c("AGEP_binned", "SEX", "RAC1P", "PUMA", "PWGTP")]
  sel_ps$S <- 0L
  
  # Align factor levels to training for selection model too
  align_levels <- function(df) {
    df$AGEP  <- factor(df$AGEP_binned,  levels = levels(nps_ca$AGEP_binned))
    df$SEX   <- factor(df$SEX,   levels = levels(nps_ca$SEX))
    df$RAC1P <- factor(df$RAC1P, levels = levels(nps_ca$RAC1P))
    # df$PUMA  <- factor(df$PUMA,  levels = PUMA_lev)
    df
  }
  sel_MR <- align_levels(sel_MR)
  sel_ps <- align_levels(sel_ps)
  
  sel_dat <- rbind(
    cbind(sel_MR, PWGTP = 1),
    sel_ps
  )
  
  if(adjust==TRUE){
    ## 1) Estimated population size from PS
    N_hat <- sum(sel_ps$PWGTP)
    ## 2) Size of nonprobability sample
    n_np <- nrow(sel_MR)
    ## 3) Adjustment factor
    adj_factor <- (N_hat - n_np) / N_hat
    ## 4) Adjust PS weights
    sel_ps$PWGTP <- sel_ps$PWGTP * adj_factor
    sel_MR$PWGTP=1
    sel_dat <- rbind(
      sel_MR[, c("AGEP_binned","SEX","RAC1P","PUMA","S","PWGTP")],
      sel_ps[, c("AGEP_binned","SEX","RAC1P","PUMA","S","PWGTP")]
    ) 
  }else{
    adj_factor=NULL
  }
  
  sel_fit <- stats::glm(
    S ~ AGEP_binned + SEX + RAC1P, 
    data = sel_dat,
    family = binomial(),
    weights = sel_dat$PWGTP
  )
  
  #bin_val is 0.11, 0.12, etc, psi is the predicted psi, using bin_val represents the random-effect group which is used during poststratification.
  predict_psi <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(list(psi = numeric(0), psi_bin = integer(0)))
    df2 <- align_levels(df)
    p <- stats::predict(sel_fit, newdata = df2, type = "response")
    p <- pmin(pmax(p, psi_eps), 1 - psi_eps)#to avoid 0
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
  
  
  
  pop_cells <- collapse_to_cells_int(
    df          = acs_pop,
    train_ref   = nps_ca,
    PUMA_lev    = PUMA_lev,
    # psi_bin_val = psi_pop$bin_val, # 0.3
    psi_vec     = psi_pop$psi,
    psi_bin_vec = psi_bin_pop, # 3
    weight_col  = "PWGTP",
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
  beta_scale <- 2.5 * sd_y / sd_x
  
  beta_scale[!is.finite(beta_scale)] <- 5
  beta_scale <- pmin(beta_scale, 5) # to avoid intercept[1]=inf
  
  beta_psi_scale <- 2.5 * sd_y / stats::sd(qlogis(psi_train$psi))
  beta_psi_scale <- min(beta_psi_scale, 5) # to avoid intercept[1]=inf
  
  sigma_psi_rate <- 1 / sd_y   
  
  
  n <- nrow(X_train)
  threads <- 4
  grainsize=as.integer(n / (threads*4))
  
  
  
  stan_data <- list(
    n = n,
    k = k,
    X = X_train,
    y = y,
    grainsize = grainsize,
    lp_psi = qlogis(psi_train$psi),                  
    G = G,
    psi_bin = psi_bin_train,
    beta_scale = as.vector(beta_scale),
    beta_psi_scale = beta_psi_scale,
    sigma_psi_rate = sigma_psi_rate
  )
  ## 3) FIT STAN (MRP-INT) ----------------------------------------------------
  fit <- mod$sample(
    data = stan_data,
    chains = 2,
    parallel_chains = 2,
    threads_per_chain = threads,
    iter_warmup = 500,
    iter_sampling = 1000,
    seed = 123
  )
  
  sum_tbl <- fit$summary()
  
  
  #using stan_glmer
  # psi_bin_train <- as.integer(bin_map[as.character(psi_train$bin_val)])
  # 
  # # logit(psi) for training
  # lp_psi_train <- qlogis(psi_train$psi)
  # 
  # ## 3) BUILD TRAINING DATA FRAME FOR stan_glmer ---------------------------
  # df_train <- nps_ca
  # df_train$PUBCOV    <- as.integer(df_train$PUBCOV)   # 0/1
  # df_train$lp_psi    <- lp_psi_train
  # df_train$psi_bin   <- factor(psi_bin_train)
  # 
  # 
  # fit_glmer <- stan_glmer(
  #   PUBCOV ~ 0 + AGEP_binned + SEX + RAC1P + lp_psi + (1 | psi_bin),
  #   data   = df_train,
  #   family = binomial(link = "logit"),
  #   chains = chains,
  #   cores  = cores,
  #   iter   = iter,
  #   warmup = warmup,
  #   seed   = seed,
  #   adapt_delta = 0.9,
  #   QR = TRUE      
  # )
  # 
  
  
  ## 4) COLLECT DRAWS & SUMMARIZE --------------------------------------------
  
  nps_ca <- MR
  nps_ca$AGEP_binned  <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX   <- factor(nps_ca$SEX)
  nps_ca$RAC1P <- factor(nps_ca$RAC1P)
  PUMA_lev     <- levels(factor(nps_ca$PUMA))
  J            <- length(PUMA_lev)
  
  params <- get_params_int_draws(fit)
  beta     <- params$beta      # S x K
  beta_psi <- params$beta_psi  # length-S
  zeta     <- params$zeta      # S x G
  
  
  
  # Pop poststrat (MRP-P with INT correction)
  mu_pop_draws <- if (nrow(pop_cells$Xp) > 0) {
    poststrat_int_SxJ(
      beta     = beta,
      beta_psi = beta_psi,
      zeta     = zeta,
      Xp       = pop_cells$Xp,
      g        = pop_cells$g,
      psi      = pop_cells$psi,
      psi_bin  = pop_cells$psi_bin,
      w        = pop_cells$w,
      J        = J
    )
  } else {
    matrix(NA_real_, nrow = nrow(beta), ncol = J)
  }
  
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-p-int")
  
  list(
    fit          = fit,
    puma_summary_mrpp = puma_summary_mrpp,
    psi_train     = psi_train,
    rhat         = sum_tbl$rhat,
    adj_factor   = adj_factor   
  )
}
