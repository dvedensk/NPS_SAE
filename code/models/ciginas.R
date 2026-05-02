ciginas_domain_mean <- function(df, response_var, weight_var) {
  design <- survey::svydesign(
    ids = ~1,
    weights = stats::as.formula(paste0("~", weight_var)),
    data = df
  )

  survey::svyby(
    stats::as.formula(paste0("~", response_var)),
    ~PUMA,
    design,
    survey::svymean,
    na.rm = TRUE,
    vartype = "se",
    keep.names = FALSE
  ) %>%
    dplyr::arrange(PUMA) %>%
    dplyr::transmute(
      PUMA,
      point_est = .data[[response_var]],
      se = se,
      variance = se^2
    )
}


ciginas_out <- function(ps,
                        nps,
                        X_formula,
                        response,
                        ps_weight_var = NULL,
                        propensity_method = "weighted") {
  if (!identical(propensity_method, "weighted")) {
    stop("The simplified Ciginas implementation currently supports only the weighted logistic IPW path.")
  }

  if (is.null(ps_weight_var)) {
    ps_weight_var <- if ("weights" %in% names(ps)) "weights" else "PWGTP"
  }

  if (!ps_weight_var %in% names(ps)) {
    stop("PS weight column not found: ", ps_weight_var)
  }

  if (!response %in% names(ps) || !response %in% names(nps)) {
    stop("Response variable must be present in both PS and NPS data frames.")
  }

  response_values <- c(ps[[response]], nps[[response]])
  bounded_response <- all(stats::na.omit(response_values) >= 0 & stats::na.omit(response_values) <= 1)

  ps_for_ipw <- ps
  ps_for_ipw$PWGTP <- ps_for_ipw[[ps_weight_var]]

  ipw_weights <- estimate_ipw(ps_for_ipw, nps, X_formula, propensity_method)
  if (any(!is.finite(ipw_weights$nps_ipw)) || any(ipw_weights$nps_ipw <= 0)) {
    stop("NPS IPW weights must be finite and strictly positive.")
  }

  ps_weighted <- ps %>% dplyr::mutate(weights = .data[[ps_weight_var]])
  nps_weighted <- nps %>% dplyr::mutate(weights = ipw_weights$nps_ipw)

  ps_est <- ciginas_domain_mean(ps_weighted, response, "weights") %>%
    dplyr::rename(
      ps_est = point_est,
      ps_se = se,
      ps_var = variance
    )

  nps_ipw_est <- ciginas_domain_mean(nps_weighted, response, "weights") %>%
    dplyr::rename(
      nps_ipw_est = point_est,
      nps_ipw_se = se,
      nps_ipw_var = variance
    )

  out <- dplyr::full_join(ps_est, nps_ipw_est, by = "PUMA") %>%
    dplyr::mutate(
      # The discrepancy term is only defined when both component estimates exist.
      nps_bias_sq = dplyr::if_else(
        is.na(ps_est) | is.na(nps_ipw_est),
        NA_real_,
        (nps_ipw_est - ps_est)^2
      ),
      nps_mse = dplyr::coalesce(nps_ipw_var, 0) + dplyr::coalesce(nps_bias_sq, 0),
      # Handle the five cases for the NPS blend weight:
      # - only NPS exists -> use NPS entirely
      # - only PS exists -> use PS entirely
      # - neither exists -> leave missing
      # - both exist but the denominator is degenerate -> fall back to PS
      # - otherwise use the usual PS variance / (PS variance + NPS MSE) blend
      weight_nps = dplyr::case_when(
        is.na(ps_est) & !is.na(nps_ipw_est) ~ 1,
        !is.na(ps_est) & is.na(nps_ipw_est) ~ 0,
        is.na(ps_est) & is.na(nps_ipw_est) ~ NA_real_,
        (ps_var + nps_mse) <= 0 ~ 0,
        TRUE ~ ps_var / (ps_var + nps_mse)
      ),
      weight_nps = pmin(pmax(weight_nps, 0), 1),
      point_est = dplyr::case_when(
        is.na(ps_est) ~ nps_ipw_est,
        is.na(nps_ipw_est) ~ ps_est,
        TRUE ~ weight_nps * nps_ipw_est + (1 - weight_nps) * ps_est
      ),
      # When both components exist, keep the blended estimate between them.
      point_est = dplyr::case_when(
        is.na(ps_est) ~ point_est,
        is.na(nps_ipw_est) ~ point_est,
        TRUE ~ pmin(
          pmax(point_est, pmin(ps_est, nps_ipw_est)),
          pmax(ps_est, nps_ipw_est)
        )
      ),
      lower_CI = NA_real_,
      upper_CI = NA_real_,
      model = "Ciginas"
    ) %>%
    dplyr::arrange(PUMA)

  if (bounded_response) {
    out <- out %>%
      dplyr::mutate(
        # Guard against tiny numerical drift outside [0, 1] for bounded outcomes.
        ps_est = pmin(pmax(ps_est, 0), 1),
        nps_ipw_est = pmin(pmax(nps_ipw_est, 0), 1),
        point_est = pmin(pmax(point_est, 0), 1)
      )
  }

  out
}
