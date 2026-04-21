
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



# ---------- 5b) getMRP_new (no data-dependent prior; beta ~ N(0,3)) ----------
# Same structure as getMRP, but DOES NOT compute beta_scale or pass it to Stan.
# Relies on si2.stan having fixed N(0,3) priors for beta.
getMRP<- function(MR,
                       ps,
                       acs_pop,
                       mod,
                       bootstrap = FALSE,
                       L = 100,
                       threads = 4,
                       n_chains = 2,
                       seed = NULL,
                       stan_iter = 1000,
                       stan_warmup = 500) {
  
  if (bootstrap && (is.null(L) || L < 1)) {
    stop("L must be specified and >= 1 when bootstrap=TRUE")
  }
  
  # TRAINING DATA ----
  nps_ca <- MR
  nps_ca$AGEP_binned <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX        <- factor(nps_ca$SEX)
  nps_ca$RAC1P      <- factor(nps_ca$RAC1P)
  PUMA_lev <- sort(unique(c(as.character(nps_ca$PUMA),
                            as.character(ps$PUMA),
                            as.character(acs_pop$PUMA))))
  
  X_train <- model.matrix(~ 0 + AGEP_binned + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train)
  n <- nrow(X_train)
  y <- as.integer(nps_ca$PUBCOV)
  
  grainsize <- as.integer(n / threads)
  puma_id <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  J <- length(PUMA_lev)
  
  stan_data <- list(
    n = n,
    k = k,
    X = X_train,
    y = y,
    J = J,
    puma_id = puma_id,
    grainsize = grainsize
  )
  
  sample_args <- list(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_chains,
    threads_per_chain = threads,
    iter_warmup = stan_warmup,
    iter_sampling = stan_iter
  )
  if (!is.null(seed)) sample_args$seed <- seed
  
  fit <- do.call(mod$sample, sample_args)
  sum_tbl <- fit$summary()

  # Collapse PS and Pop to cells
  ps_cells  <- collapse_to_cells(ps,      nps_ca, PUMA_lev, weight_col = "weights")
  pop_cells <- collapse_to_cells(acs_pop, nps_ca, PUMA_lev, weight_col = "weights")
  
  beta <- get_beta_draws(fit)
  apuma_draws <- get_apuma_draws(fit)
  
  mu_ps_draws <- if (nrow(ps_cells$Xp) > 0) {
    poststrat_SxJ(
      beta = beta,
      a_puma = apuma_draws,
      Xp = ps_cells$Xp,
      g = ps_cells$g,
      w = ps_cells$w,
      J = J
    )
  } else {
    matrix(NA_real_, nrow = nrow(beta), ncol = J)
  }
  
  mu_pop_draws <- if (nrow(pop_cells$Xp) > 0) {
    poststrat_SxJ(
      beta = beta,
      a_puma = apuma_draws,
      Xp = pop_cells$Xp,
      g = pop_cells$g,
      w = pop_cells$w,
      J = J
    )
  } else {
    matrix(NA_real_, nrow = nrow(beta), ncol = J)
  }
  
  puma_summary_mrpr <- make_summary(mu_ps_draws,  PUMA_lev, "mrp-r")
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev, "mrp-p")
  
  puma_summary_mrpr_bootstrap <- NULL
  if (bootstrap) {
    S <- nrow(beta)
    mu_boot_array <- array(NA_real_, dim = c(S, J, L))
    cat("\nGenerating", L, "bootstrap samples for getMRP_new...\n")
    for (l in 1:L) {
      boot_idx <- sample(
        seq_len(nrow(ps)),
        size = nrow(ps),
        replace = TRUE,
        prob = ps$weights
      )
      ps_boot <- ps[boot_idx, ]
      ps_boot_cells <- collapse_to_cells(ps_boot, nps_ca, PUMA_lev, weight_col = "weights")
      mu_boot_array[, , l] <- if (nrow(ps_boot_cells$Xp) > 0) {
        poststrat_SxJ(
          beta = beta,
          a_puma = apuma_draws,
          Xp = ps_boot_cells$Xp,
          g = ps_boot_cells$g,
          w = ps_boot_cells$w,
          J = J
        )
      } else {
        matrix(NA_real_, nrow = S, ncol = J)
      }
      if (l %% 10 == 0) {
        cat("  Completed", l, "of", L, "bootstrap samples\n")
      }
    }
    mu_combined <- matrix(aperm(mu_boot_array, c(3, 1, 2)), nrow = S * L, ncol = J)
    puma_summary_mrpr_bootstrap <- make_summary(
      mu_combined, PUMA_lev, "mrp-r-bootstrap"
    )
  }
  
  list(
    puma_summary_mrpr = puma_summary_mrpr,
    puma_summary_mrpp = puma_summary_mrpp,
    puma_summary_mrpr_bootstrap = puma_summary_mrpr_bootstrap,
    rhat = sum_tbl$rhat
  )
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
                                  bin_digits = 2,
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
  cells$psi_bin <- map_bins(bin_fun(cells$psi, digits = bin_digits))
  
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


# ---------- 3b) Cell-level inclusion probability (external step) ----------
# Returns one psi per unique (AGEP_binned, SEX, RAC1P, PUMA). Optional cell_psi.
# include_response = TRUE adds PUBCOV to the selection formula; population psi is
# then marginalized over PUBCOV proportions from acs_pop.
compute_cell_inclusion_probs <- function(ps,
                                         nps,
                                         train_ref,
                                         PUMA_lev,
                                         acs_pop = NULL,
                                         cell_psi = NULL,
                                         psi_eps = 1e-6,
                                         bin_digits = 2,
                                         adjust = FALSE,
                                         include_response = FALSE) {
  base_cols <- c("AGEP_binned", "SEX", "RAC1P", "PUMA")
  sel_cols  <- if (include_response) c(base_cols, "PUBCOV") else base_cols
  align_lev <- function(df) {
    df$AGEP_binned <- factor(df$AGEP_binned, levels = levels(train_ref$AGEP_binned))
    df$SEX         <- factor(df$SEX,         levels = levels(train_ref$SEX))
    df$RAC1P       <- factor(df$RAC1P,       levels = levels(train_ref$RAC1P))
    df$PUMA        <- factor(df$PUMA,        levels = PUMA_lev)
    if (include_response && "PUBCOV" %in% names(df))
      df$PUBCOV <- factor(as.character(df$PUBCOV), levels = c("0", "1"))
    df
  }
  if (!is.null(cell_psi)) {
    cells <- align_lev(cell_psi[, c(base_cols, "psi")])
    cells$psi <- pmin(pmax(cells$psi, psi_eps), 1 - psi_eps)
    cells$bin_val <- bin_fun(cells$psi, digits = bin_digits)
    return(list(cells_nps = cells, cells_ps = cells, cells_pop = cells))
  }
  if (include_response && (!"PUBCOV" %in% names(train_ref) || !"PUBCOV" %in% names(ps)))
    stop("include_response = TRUE requires PUBCOV in both NPS and PS data.")
  sel_MR <- train_ref[, sel_cols]
  sel_MR$S <- 1L
  sel_ps  <- ps[, c(sel_cols, "weights")]
  sel_ps$S <- 0L
  sel_MR <- align_lev(sel_MR)
  sel_ps <- align_lev(sel_ps)
  if (adjust) {
    N_hat <- sum(sel_ps$weights)
    n_np  <- nrow(sel_MR)
    sel_ps$weights <- sel_ps$weights * (N_hat - n_np) / N_hat
  }
  sel_dat <- dplyr::bind_rows(
    dplyr::mutate(sel_MR[, c(sel_cols, "S")], weights = 1),
    sel_ps[, c(sel_cols, "S", "weights")]
  )
  sel_formula <- if (include_response)
    S ~ AGEP_binned + SEX + RAC1P + PUBCOV
  else
    S ~ AGEP_binned + SEX + RAC1P + PUMA
  sel_fit <- stats::glm(sel_formula, data = sel_dat, family = binomial(),
                        weights = sel_dat$weights)
  if (include_response) {
    # Per-individual psi averaged to cell level (marginalizes over observed PUBCOV)
    nps_unit <- align_lev(train_ref[, sel_cols])
    nps_unit$psi_unit <- pmin(pmax(
      stats::predict(sel_fit, newdata = nps_unit, type = "response"), psi_eps), 1 - psi_eps)
    cells_nps <- nps_unit %>%
      dplyr::group_by(AGEP_binned, SEX, RAC1P, PUMA) %>%
      dplyr::summarise(psi = mean(psi_unit, na.rm = TRUE), .groups = "drop")
    ps_unit <- align_lev(ps[, sel_cols])
    ps_unit$psi_unit <- pmin(pmax(
      stats::predict(sel_fit, newdata = ps_unit, type = "response"), psi_eps), 1 - psi_eps)
    cells_ps <- ps_unit %>%
      dplyr::group_by(AGEP_binned, SEX, RAC1P, PUMA) %>%
      dplyr::summarise(psi = mean(psi_unit, na.rm = TRUE), .groups = "drop")
    cells_pop <- NULL
    if (!is.null(acs_pop) && nrow(acs_pop) > 0 && "PUBCOV" %in% names(acs_pop)) {
      acs_lev <- .align_to_train(acs_pop, train_ref, PUMA_lev)
      acs_lev$PUBCOV <- factor(as.character(acs_lev$PUBCOV), levels = c("0", "1"))
      cell_prop <- acs_lev %>%
        dplyr::group_by(AGEP_binned, SEX, RAC1P, PUMA) %>%
        dplyr::summarise(n = dplyr::n(), n1 = sum(PUBCOV == "1", na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(prop1 = n1 / n)
      cells_grid <- align_lev(cell_prop[, base_cols])
      cells_grid$PUBCOV <- factor("0", levels = c("0", "1"))
      psi0 <- pmin(pmax(stats::predict(sel_fit, newdata = cells_grid, type = "response"), psi_eps), 1 - psi_eps)
      cells_grid$PUBCOV <- factor("1", levels = c("0", "1"))
      psi1 <- pmin(pmax(stats::predict(sel_fit, newdata = cells_grid, type = "response"), psi_eps), 1 - psi_eps)
      cells_pop <- cell_prop %>%
        dplyr::mutate(psi = as.numeric((1 - prop1) * psi0 + prop1 * psi1))
      cells_pop <- align_lev(cells_pop[, c(base_cols, "psi")])
    }
  } else {
    cells_nps <- dplyr::distinct(train_ref[, base_cols])
    cells_ps  <- dplyr::distinct(ps[, base_cols])
    cells     <- dplyr::distinct(dplyr::bind_rows(cells_nps, cells_ps))
    if (!is.null(acs_pop) && nrow(acs_pop) > 0) {
      acs_lev <- .align_to_train(acs_pop, train_ref, PUMA_lev)
      cells_acs <- dplyr::distinct(acs_lev[, base_cols])
      cells <- dplyr::distinct(dplyr::bind_rows(cells, cells_acs))
    }
    cells <- align_lev(cells)
    p <- stats::predict(sel_fit, newdata = cells, type = "response")
    p <- pmin(pmax(p, psi_eps), 1 - psi_eps)
    cells$psi <- as.numeric(p)
    cells_nps <- cells
    cells_ps  <- cells
    cells_pop <- cells
  }
  cells_nps$bin_val <- bin_fun(cells_nps$psi, digits = bin_digits)
  cells_ps$bin_val  <- bin_fun(cells_ps$psi,  digits = bin_digits)
  if (!is.null(cells_pop)) cells_pop$bin_val <- bin_fun(cells_pop$psi, digits = bin_digits)
  list(cells_nps = cells_nps, cells_ps = cells_ps, cells_pop = cells_pop)
}




# ---------- 4b) getMRP_INT_new (cell-level psi ) ----------
getMRP_INT <- function(MR,
                       ps,
                       acs_pop,
                       mod,
                       cell_psi = NULL,
                       include_response = FALSE,
                       bootstrap = FALSE,
                       L = 100,
                       threads = 4,
                       n_chains = 2,
                       seed = NULL,
                       stan_iter = 1000,
                       stan_warmup = 500,
                       adjust = FALSE,
                       psi_eps = 1e-6,
                       bin_digits = 2) {
  if (bootstrap && (is.null(L) || L < 1)) stop("L must be specified and >= 1 when bootstrap=TRUE")
  nps_ca <- MR
  nps_ca$AGEP_binned <- factor(nps_ca$AGEP_binned)
  nps_ca$SEX        <- factor(nps_ca$SEX)
  nps_ca$RAC1P      <- factor(nps_ca$RAC1P)
  PUMA_lev <- sort(unique(c(
    as.character(nps_ca$PUMA),
    as.character(ps$PUMA),
    as.character(acs_pop$PUMA)
  )))
  J <- length(PUMA_lev)
  X_train <- model.matrix(~ 0 + AGEP_binned + SEX + RAC1P, data = nps_ca)
  k <- ncol(X_train)
  n <- nrow(X_train)
  y <- as.integer(nps_ca$PUBCOV)
  train_ref <- nps_ca
  cells_psi <- compute_cell_inclusion_probs(
    ps = ps, nps = nps_ca, train_ref = train_ref, PUMA_lev = PUMA_lev,
    acs_pop = acs_pop, cell_psi = cell_psi, psi_eps = psi_eps, bin_digits = bin_digits,
    adjust = adjust, include_response = include_response
  )
  cells_nps <- cells_psi$cells_nps
  cells_ps  <- cells_psi$cells_ps
  cells_pop <- cells_psi$cells_pop
  if (is.null(cells_pop) && !is.null(acs_pop) && nrow(acs_pop) > 0) {
    if (isTRUE(include_response)) {
      stop(
        "include_response = TRUE requires acs_pop with PUBCOV so population ",
        "PUBCOV composition can be marginalized."
      )
    }
    cells_pop <- cells_nps
  }
  all_bin_vals <- sort(unique(c(
    cells_nps$bin_val,
    cells_ps$bin_val,
    if (!is.null(cells_pop)) cells_pop$bin_val else NULL
  )))
  G <- length(all_bin_vals)
  bin_map <- setNames(seq_along(all_bin_vals), as.character(all_bin_vals))
  map_bins <- function(bin_val_vec) as.integer(bin_map[as.character(bin_val_vec)])
  cells_nps$psi_bin <- map_bins(cells_nps$bin_val)
  cells_ps$psi_bin  <- map_bins(cells_ps$bin_val)
  if (!is.null(cells_pop)) cells_pop$psi_bin <- map_bins(cells_pop$bin_val)
  nps_ca$PUMA <- factor(nps_ca$PUMA, levels = PUMA_lev)
  nps_join <- dplyr::left_join(
    nps_ca[, c("AGEP_binned", "SEX", "RAC1P", "PUMA")],
    cells_nps[, c("AGEP_binned", "SEX", "RAC1P", "PUMA", "psi", "psi_bin")],
    by = c("AGEP_binned", "SEX", "RAC1P", "PUMA")
  )
  if (any(is.na(nps_join$psi))) stop("Some NPS cells missing in cell_psi.")
  lp_psi_train  <- qlogis(pmin(pmax(nps_join$psi, psi_eps), 1 - psi_eps))
  psi_bin_train <- nps_join$psi_bin
  ps_align <- ps
  ps_align$AGEP_binned <- factor(ps_align$AGEP_binned, levels = levels(train_ref$AGEP_binned))
  ps_align$SEX         <- factor(ps_align$SEX,         levels = levels(train_ref$SEX))
  ps_align$RAC1P      <- factor(ps_align$RAC1P,      levels = levels(train_ref$RAC1P))
  ps_align$PUMA       <- factor(ps_align$PUMA,      levels = PUMA_lev)
  ps_join <- dplyr::left_join(
    ps_align[, c("AGEP_binned", "SEX", "RAC1P", "PUMA")],
    cells_ps[, c("AGEP_binned", "SEX", "RAC1P", "PUMA", "psi", "psi_bin")],
    by = c("AGEP_binned", "SEX", "RAC1P", "PUMA")
  )
  psi_ps     <- if (any(is.na(ps_join$psi))) rep(NA_real_, nrow(ps)) else ps_join$psi
  psi_bin_ps <- if (any(is.na(ps_join$psi_bin))) rep(NA_integer_, nrow(ps)) else ps_join$psi_bin
  acs_align <- .align_to_train(acs_pop, train_ref, PUMA_lev)
  cells_pop_use <- if (is.null(cells_pop)) cells_nps else cells_pop
  acs_join <- dplyr::left_join(
    acs_align[, c("AGEP_binned", "SEX", "RAC1P", "PUMA")],
    cells_pop_use[, c("AGEP_binned", "SEX", "RAC1P", "PUMA", "psi", "psi_bin")],
    by = c("AGEP_binned", "SEX", "RAC1P", "PUMA")
  )
  psi_pop     <- acs_join$psi
  psi_bin_pop <- acs_join$psi_bin
  if (any(is.na(psi_pop))) {
    fill_psi <- mean(cells_nps$psi, na.rm = TRUE)
    psi_pop[is.na(psi_pop)]     <- fill_psi
    psi_bin_pop[is.na(psi_bin_pop)] <- map_bins(bin_fun(fill_psi, digits = bin_digits))
  }
  pop_cells <- collapse_to_cells_int(
    df = acs_pop, train_ref = train_ref, PUMA_lev = PUMA_lev,
    psi_vec = psi_pop, psi_bin_vec = psi_bin_pop,
    weight_col = "weights", bin_map = bin_map, bin_digits = bin_digits, map_bins = map_bins
  )
  ps_cells <- collapse_to_cells_int(
    df = ps, train_ref = train_ref, PUMA_lev = PUMA_lev,
    psi_vec = psi_ps, psi_bin_vec = psi_bin_ps,
    weight_col = "weights", bin_map = bin_map, bin_digits = bin_digits, map_bins = map_bins
  )
  grainsize <- max(1L, as.integer(n / (threads * 4)))
  puma_id_train <- as.integer(factor(nps_ca$PUMA, levels = PUMA_lev))
  sd_y <- stats::sd(as.numeric(y))
  if (!is.finite(sd_y) || sd_y <= 0) sd_y <- 1
  stan_data <- list(
    n = n, k = k, X = X_train, y = y,
    grainsize = grainsize,
    lp_psi = as.vector(lp_psi_train),
    G = G, psi_bin = as.integer(psi_bin_train),
    J = J, puma_id = puma_id_train,
    sigma_psi_rate = 1 / sd_y
  )
  sample_args <- list(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_chains,
    threads_per_chain = threads,
    iter_warmup = stan_warmup,
    iter_sampling = stan_iter
  )
  if (!is.null(seed)) sample_args$seed <- seed
  fit <- do.call(mod$sample, sample_args)
  sum_tbl <- fit$summary()
  params <- get_params_int_draws(fit)
  beta <- params$beta
  beta_psi <- params$beta_psi
  zeta <- params$zeta
  a_puma <- params$a_puma
  mu_pop_draws <- if (nrow(pop_cells$Xp) > 0) {
    poststrat_int_SxJ(beta, beta_psi, zeta, a_puma,
                      pop_cells$Xp, pop_cells$g, pop_cells$psi, pop_cells$psi_bin,
                      pop_cells$w, J)
  } else matrix(NA_real_, nrow = nrow(beta), ncol = J)
  mu_ps_draws <- if (nrow(ps_cells$Xp) > 0) {
    poststrat_int_SxJ(beta, beta_psi, zeta, a_puma,
                      ps_cells$Xp, ps_cells$g, ps_cells$psi, ps_cells$psi_bin,
                      ps_cells$w, J)
  } else matrix(NA_real_, nrow = nrow(beta), ncol = J)
  lbl <- if (include_response) "-PUBCOV" else ""
  pumas_pop <- sort(unique(pop_cells$g))
  pumas_ps  <- sort(unique(ps_cells$g))
  PUMA_lev_mrpp <- if (length(pumas_pop) < J) PUMA_lev[pumas_pop] else PUMA_lev
  PUMA_lev_mrpr <- if (length(pumas_ps)  < J) PUMA_lev[pumas_ps]  else PUMA_lev
  if (length(pumas_pop) < J) mu_pop_draws <- mu_pop_draws[, pumas_pop, drop = FALSE]
  if (length(pumas_ps)  < J) mu_ps_draws  <- mu_ps_draws[,  pumas_ps,  drop = FALSE]
  puma_summary_mrpp <- make_summary(mu_pop_draws, PUMA_lev_mrpp, paste0("mrp-p-int", lbl))
  puma_summary_mrpr <- make_summary(mu_ps_draws,  PUMA_lev_mrpr, paste0("mrp-r-int", lbl))
  puma_summary_mrpr_bootstrap <- NULL
  if (bootstrap) {
    S <- nrow(beta)
    mu_boot_array <- array(NA_real_, dim = c(S, J, L))
    cat("\nGenerating", L, "bootstrap samples for getMRP_INT...\n")
    for (l in 1:L) {
      boot_idx <- sample(seq_len(nrow(ps)), nrow(ps), replace = TRUE, prob = ps$weights)
      ps_boot <- ps[boot_idx, ]
      ps_boot$AGEP_binned <- factor(ps_boot$AGEP_binned, levels = levels(train_ref$AGEP_binned))
      ps_boot$SEX         <- factor(ps_boot$SEX,         levels = levels(train_ref$SEX))
      ps_boot$RAC1P      <- factor(ps_boot$RAC1P,      levels = levels(train_ref$RAC1P))
      ps_boot$PUMA       <- factor(ps_boot$PUMA,       levels = PUMA_lev)
      ps_boot_join <- dplyr::left_join(
        ps_boot[, c("AGEP_binned", "SEX", "RAC1P", "PUMA")],
        cells_ps[, c("AGEP_binned", "SEX", "RAC1P", "PUMA", "psi", "psi_bin")],
        by = c("AGEP_binned", "SEX", "RAC1P", "PUMA")
      )
      psi_boot <- ps_boot_join$psi
      psi_bin_boot <- ps_boot_join$psi_bin
      if (any(is.na(psi_boot))) psi_boot[is.na(psi_boot)] <- mean(cells_nps$psi)
      if (any(is.na(psi_bin_boot))) psi_bin_boot[is.na(psi_bin_boot)] <- 1L
      ps_boot_cells <- collapse_to_cells_int(
        df = ps_boot, train_ref = train_ref, PUMA_lev = PUMA_lev,
        psi_vec = psi_boot, psi_bin_vec = psi_bin_boot,
        weight_col = "weights", bin_map = bin_map, bin_digits = bin_digits, map_bins = map_bins
      )
      mu_boot_array[, , l] <- if (nrow(ps_boot_cells$Xp) > 0) {
        poststrat_int_SxJ(beta, beta_psi, zeta, a_puma,
                          ps_boot_cells$Xp, ps_boot_cells$g, ps_boot_cells$psi, ps_boot_cells$psi_bin,
                          ps_boot_cells$w, J)
      } else matrix(NA_real_, S, J)
      if (l %% 10 == 0) cat("  Completed", l, "of", L, "bootstrap samples\n")
    }
    mu_combined <- matrix(aperm(mu_boot_array, c(3, 1, 2)), nrow = S * L, ncol = J)
    if (length(pumas_ps) < J) mu_combined <- mu_combined[, pumas_ps, drop = FALSE]
    puma_summary_mrpr_bootstrap <- make_summary(mu_combined, PUMA_lev_mrpr, paste0("mrp-r-int", lbl, "-bootstrap"))
  }
  list(
    fit = fit,
    puma_summary_mrpp = puma_summary_mrpp,
    puma_summary_mrpr = puma_summary_mrpr,
    puma_summary_mrpr_bootstrap = puma_summary_mrpr_bootstrap,
    rhat = sum_tbl$rhat,
    cell_psi = cells_nps
  )
}
