## helpers for outcome model
make_mrp_inclusion_cells <- function(pop_df, nps, cell_vars) {
  pop_cells <- pop_df %>%
    filter(!if_any(all_of(cell_vars), is.na)) %>%
    count(across(all_of(cell_vars)), name = "N")

  nps_cells <- nps %>%
    filter(!if_any(all_of(cell_vars), is.na)) %>%
    count(across(all_of(cell_vars)), name = "n")

  pop_cells %>%
    left_join(nps_cells, by = cell_vars) %>%
    mutate(n = replace_na(n, 0L)) %>%
    arrange(across(all_of(cell_vars)))
}

make_mrp_inclusion_stan_data <- function(cells, X_formula, puma_random=FALSE) {
  #drop intercept if it exists
  X <- model.matrix(X_formula, data = cells)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  stan_data <- list(
    J = nrow(cells),
    n = as.integer(cells$n),
    N = as.numeric(cells$N),
    K = ncol(X),
    X = X
  )

  if (puma_random) {
    puma_id <- as.integer(factor(cells$PUMA))
    stan_data$G_puma <- max(puma_id)
    stan_data$puma_id <- puma_id
  }

  return(stan_data)  
}

get_ipw_weights <- function(nps, cells, cell_vars) {
  joined <- nps %>%
    mutate(.row_id = row_number()) %>%
    left_join(
      cells %>%
        select(all_of(cell_vars), ipw),
      by = cell_vars
    ) %>%
    arrange(.row_id)

  return(joined$ipw)
}

fit_mrp_inclusion_helper <- function(cells, X_formula, mod,
                                     puma_random = FALSE,
                                     nps = NULL,
                                     cell_vars = NULL,
                                     response_var = NULL,
                                     model_name = NULL,
                                     chains = 4, iter_warmup = 1000,
                                     iter_sampling = 1000, seed = 123,
                                     threads_per_chain = 1) {
  stan_data <- make_mrp_inclusion_stan_data(
    cells = cells,
    X_formula = X_formula,
    puma_random = puma_random
  )

  fit <- mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = 250
  )

  psi_draws <- as_draws_matrix(fit$draws("psi"))
  logit_psi_draws <- as_draws_matrix(fit$draws("logit_psi"))

  cells_out <- cells %>%
    mutate(
      psi_hat = apply(psi_draws, 2, median),
      logit_psi_hat = apply(logit_psi_draws, 2, median),
      psi_group = factor(round(psi_hat, 1)),
      ipw = 1 / psi_hat
    )

  weights <- NULL
  ipw_HT <- NULL

  if (!is.null(nps) && !is.null(cell_vars)) {
    weights <- get_ipw_weights(
      nps = nps,
      cells = cells_out,
      cell_vars = cell_vars
    )
  }

  if (!is.null(weights) && !is.null(response_var) && !is.null(model_name)) {
    ipw_df <- data.frame(
      response = nps[[response_var]],
      PUMA = nps$PUMA,
      weights = weights
    )

    names(ipw_df)[1] <- response_var

    ipw_HT <- HT(ipw_df, model_name, response_var = response_var)
  }

  return(list(
    fit = fit,
    cells = cells_out,
    stan_data = stan_data,
    weights = weights,
    ipw_HT = ipw_HT
  ))
}

fit_mrp_inclusion <- function(cells, X_formula,
                              puma_mode = c("none", "fixed", "random"),
                              fixed_mod, random_mod,
                              nps = NULL,
                              cell_vars = NULL,
                              response_var = NULL,
                              model_name = NULL,
                              chains = 2, iter_warmup = 250,
                              iter_sampling = 250, seed = 123,
                              threads_per_chain = 1) {
  puma_mode <- match.arg(puma_mode)

  if (puma_mode == "none") {
    return(
      fit_mrp_inclusion_helper(
        cells = cells,
        X_formula = X_formula,
        mod = fixed_mod,
        puma_random = FALSE,
        nps = nps,
        cell_vars = cell_vars,
        response_var = response_var,
        model_name = model_name,
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        seed = seed,
        threads_per_chain = threads_per_chain
      )
    )
  }

  if (puma_mode == "fixed") {
    return(
      fit_mrp_inclusion_helper(
        cells = cells,
        X_formula = update(X_formula, ~ . + PUMA),
        mod = fixed_mod,
        puma_random = FALSE,
        nps = nps,
        cell_vars = cell_vars,
        response_var = response_var,
        model_name = model_name,
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        seed = seed,
        threads_per_chain = threads_per_chain
      )
    )
  }

  fit_mrp_inclusion_helper(
    cells = cells,
    X_formula = X_formula,
    mod = random_mod,
    puma_random = TRUE,
    nps = nps,
    cell_vars = cell_vars,
    response_var = response_var,
    model_name = model_name,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    threads_per_chain = threads_per_chain
  )
}

## helpers for outcome model
make_mrp_outcome_cells <- function(cells, nps, response, cell_vars) {
  nps_y_cells <- nps %>%
    filter(!if_any(all_of(c(cell_vars, response)), is.na)) %>%
    group_by(across(all_of(cell_vars))) %>%
    summarise(
      m = n(),
      y = sum(.data[[response]]),
      .groups = "drop"
    )

  cells %>%
    left_join(nps_y_cells, by = cell_vars) %>%
    mutate(
      m = replace_na(m, 0L),
      y = replace_na(y, 0L),
      failures = m - y
    )
}

make_mrp_outcome_stan_data <- function(cells, X_formula) {
  X <- model.matrix(X_formula, data = cells)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  puma_factor <- factor(cells$PUMA)
  psi_group_factor <- factor(cells$psi_group)

  puma_id <- as.integer(puma_factor)
  psi_group_id <- as.integer(psi_group_factor)

  G_puma <- max(puma_id)
  J <- nrow(cells)

  A_puma <- matrix(0, nrow = G_puma, ncol = J)
  #place 1s connecting puma_id to cell number
  A_puma[cbind(puma_id, seq_len(J))] <- 1 

  obs_id <- which(cells$m > 0)

  stan_data <- list(
    J = J,
    J_obs = length(obs_id),
    obs_id = as.integer(obs_id),
    y = as.integer(cells$y[obs_id]),
    m = as.integer(cells$m[obs_id]),
    K = ncol(X),
    X = X,
    logit_psi = as.numeric(cells$logit_psi_hat),
    G_puma = G_puma,
    puma_id = puma_id,
    G_psi = max(psi_group_id),
    psi_group_id = psi_group_id,
    N = as.numeric(cells$N),
    A_puma = A_puma
  )

  list(
    stan_data = stan_data,
    puma_levels = levels(puma_factor),
    X = X
  )
}

fit_mrp_outcome <- function(cells, X_formula, mod,
                            chains = 2, iter_warmup = 250,
                            iter_sampling = 250, seed = 123,
                            threads_per_chain = 1) {
  data_obj <- make_mrp_outcome_stan_data(
    cells = cells,
    X_formula = X_formula
  )

  fit <- mod$sample(
    data = data_obj$stan_data,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = 250
  )

  ybar_puma_draws <- as_draws_matrix(fit$draws("ybar_puma"))

  summary_df <- tibble(
    PUMA = data_obj$puma_levels,
    point_est = apply(ybar_puma_draws, 2, mean),   
    lower_CI = apply(ybar_puma_draws, 2, quantile, probs = 0.025),
    upper_CI = apply(ybar_puma_draws, 2, quantile, probs = 0.975)
  )

  return(list(
    fit = fit,
    cells = cells,
    stan_data = data_obj$stan_data,
    summary_df = summary_df
  ))
}
