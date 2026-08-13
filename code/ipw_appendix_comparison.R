library(doParallel)
library(dplyr)
library(foreach)
library(ggplot2)
library(knitr)
library(kableExtra)
library(lme4)
library(nonprobsampling)
library(readr)
library(sampling)
library(stringr)
library(survey)
library(tidyr)
library(tidyverse)

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R"))

response_var <- "PUBCOV"
response_type <- "binary"

settings <- c("easy", "medium", "hard")
Nsim <- 100
alpha <- 0.05
numCores <- 10

PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)

get_nps_config <- function(setting) {
  if (setting == "easy") {
    list(PWGTP = 0.10, POVPIP = -0.22)
  } else if (setting == "medium") {
    list(PWGTP = 0.10, POVPIP = -0.41)
  } else if (setting == "hard") {
    list(PWGTP = 0.10, POVPIP = -1.52)
  } else {
    stop("Improper setting")
  }
}

get_save_dir <- function(setting) {
  file.path("ipw_data", paste0(setting, "_setting_pubcov_ipw_test"))
}

acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv"), show_col_types = FALSE) %>%
  mutate(
    AGEP_binned = factor(AGEP_binned),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    response_true = mean(.data[[response_var]], na.rm = TRUE),
    .groups = "drop"
  )

formulas_to_run <- list(
  agep = as.formula("~ AGEP + RAC1P + SEX"),
  agep_y = as.formula(paste0("~ AGEP + RAC1P + SEX + ", response_var)),
  agep_binned = as.formula("~ AGEP_binned + RAC1P + SEX"),
  agep_binned_y = as.formula(paste0("~ AGEP_binned + RAC1P + SEX + ", response_var))
)

formula_labels <- c(
  agep = "AGEP",
  agep_y = "AGEP+Y",
  agep_binned = "AGEP_binned",
  agep_binned_y = "AGEP_binned+Y"
)

ev_formula <- as.formula("~ AGEP_binned + RAC1P + SEX")
ev_formula_y <- as.formula(paste0("~ AGEP_binned + RAC1P + SEX + ", response_var))

nonprobsampling_ipw_wrapper <- function(method,
                                        formula_name,
                                        nps,
                                        samp_design,
                                        X_formula,
                                        response_var) {
  method_labels <- c(
    alp = "IPW (ALP)",
    clw = "IPW (CLW)",
    calibration = "IPW (calibration)"
  )

  fit <- est_pw(
    data = list(nps, samp_design),
    p_formula = X_formula,
    method = method,
    control = pw_solver_control(ftol = 1e-6)
  )

  out <- pwmean(fit, y = response_var, zcol = "PUMA")

  label <- paste0(
    method_labels[[method]],
    " [",
    formula_labels[[formula_name]],
    "]"
  )

  domain_est <- out$estimates %>%
    transmute(
      PUMA = str_remove(domain, "^PUMA =\\s*"),
      point_est = adjusted_mean,
      lower_CI = adjusted_lower,
      upper_CI = adjusted_upper,
      model = label
    )

  list(
    method = method,
    formula_name = formula_name,
    label = label,
    fit = fit,
    pwmean = out,
    weights = fit$pseudo_weights,
    domain_est = domain_est
  )
}

run_ev_estimator <- function(ps,
                             nps,
                             y,
                             PUMAs,
                             formula,
                             method,
                             model_label,
                             response_var) {
  weights_out <- estimate_ipw(ps, nps, formula, method)
  weights <- c(weights_out$ps_ipw, weights_out$nps_ipw)

  df <- data.frame(
    response = y,
    PUMA = PUMAs,
    weights = weights
  )

  names(df)[1] <- response_var

  out <- HT(df, model_label, response_var = response_var)
  out$model <- model_label

  out
}

all_results <- list()
all_summaries <- list()

for (setting in settings) {
  save_dir <- get_save_dir(setting)
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  NPS_weight_config <- get_nps_config(setting)
  N_pop <- nrow(acs_pop)
  pop_mean <- mean(acs_pop[[response_var]], na.rm = TRUE)

  cl <- makeCluster(numCores, type = "FORK", outfile = "")
  registerDoParallel(cl)

  setting_results <- foreach(
    sim = seq_len(Nsim),
    .combine = "rbind",
    .verbose = TRUE,
    .packages = c("tidyverse", "survey", "nonprobsampling", "lme4", "betareg")
  ) %dopar% {
    curr_seed <- 99 + sim
    set.seed(curr_seed)

    print(paste0("Starting setting", setting, "simulation", sim))

    ps_sample <- get_strat_PS(
      pop_df = acs_pop,
      samp_frac = 0.005,
      weight_config = PS_weight_config
    )

    nps_sample <- get_NPS(
      pop_df = acs_pop,
      samp_frac = 0.05,
      weight_config = NPS_weight_config,
      internet_only = FALSE
    )

    ps <- acs_pop[ps_sample$idx, ] %>%
      mutate(weights = ps_sample$weights)

    nps <- acs_pop[nps_sample$idx, ]

    samp_design <- svydesign(
      ids = ~1,
      strata = ~PUMA,
      weights = ~weights,
      data = ps
    )

    direst <- svyby(
      as.formula(paste0("~", response_var)),
      ~PUMA,
      samp_design,
      svymean,
      na.rm = TRUE,
      vartype = "se",
      keep.names = FALSE
    ) %>%
      arrange(PUMA) %>%
      transmute(
        PUMA,
        point_est = .data[[response_var]],
        lower_CI = point_est + qnorm(alpha / 2) * se,
        upper_CI = point_est + qnorm(1 - alpha / 2) * se,
        model = "direst"
      )

    ipw_methods <- c("alp", "clw", "calibration")
    ipw_results <- list()

    for (formula_name in names(formulas_to_run)) {
      for (method in ipw_methods) {
        key <- paste(method, formula_name, sep = "_")

        ipw_results[[key]] <- nonprobsampling_ipw_wrapper(
          method = method,
          formula_name = formula_name,
          nps = nps,
          samp_design = samp_design,
          X_formula = formulas_to_run[[formula_name]],
          response_var = response_var
        )
      }
    }

    ipw_domain_est <- bind_rows(lapply(ipw_results, `[[`, "domain_est"))

    y_ps <- ps[[response_var]]
    y_nps <- nps[[response_var]]
    y <- c(y_ps, y_nps)
    PUMAs <- c(ps$PUMA, nps$PUMA)

    ev_HT <- run_ev_estimator(
      ps = ps,
      nps = nps,
      y = y,
      PUMAs = PUMAs,
      formula = ev_formula,
      method = "weighted",
      model_label = "IPW (E&V)",
      response_var = response_var
    )

    ev_y_HT <- run_ev_estimator(
        ps = ps,
        nps = nps,
        y = y,
        PUMAs = PUMAs,
        formula = ev_formula_y,
        method = "weighted",
        model_label = "IPW (E&V)+Y",
        response_var = response_var
    )

    ev_re_HT <- run_ev_estimator(
        ps = ps,
        nps = nps,
        y = y,
        PUMAs = PUMAs,
        formula = ev_formula,
        method = "reff",
        model_label = "IPW (E&V RE)",
        response_var = response_var
    )

    ev_re_y_HT <- run_ev_estimator(
        ps = ps,
        nps = nps,
        y = y,
        PUMAs = PUMAs,
        formula = ev_formula_y,
        method = "reff",
        model_label = "IPW (E&V RE)+Y",
        response_var = response_var
    )

    results <- bind_rows(
        direst,
        ipw_domain_est,
        ev_HT,
        ev_y_HT,
        ev_re_HT,
        ev_re_y_HT
    )

    results$sim_num <- sim
    results$setting <- setting

    save(
      results,
      file = file.path(save_dir, paste0("intermediate_result_", sim, ".RData"))
    )

    results
  }

  stopCluster(cl)
  registerDoSEQ()

  setting_results <- tibble(setting_results)

  summary_df <- true_values %>%
    left_join(setting_results, by = "PUMA") %>%
    group_by(setting, model) %>%
    summarize(
      MSE = mean((response_true - point_est)^2, na.rm = TRUE),
      MAB = mean(abs(response_true - point_est), na.rm = TRUE),
      Coverage = mean(between(response_true, lower_CI, upper_CI), na.rm = TRUE),
      `Int. Score` = mean(int_score(alpha, response_true, lower_CI, upper_CI), na.rm = TRUE),
      .groups = "drop"
    )

  save(
    setting_results,
    summary_df,
    file = file.path(save_dir, "ACS_NPS_ipw_test_results.RData")
  )

  print(summary_df %>% arrange(MSE))

  all_results[[setting]] <- setting_results
  all_summaries[[setting]] <- summary_df
}

results_df <- bind_rows(all_results)
summary_df <- bind_rows(all_summaries)

save(
  results_df,
  summary_df,
  file = file.path("ipw_data", "ACS_NPS_ipw_test_all_settings.RData")
)

print(summary_df %>% arrange(setting, MSE))


### figures
alpha <- 0.05
response_var <- "PUBCOV"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

setting_dirs <- tibble(
  setting_raw = c("easy", "medium", "hard"),
  setting = c("Favorable", "Typical", "Extreme"),
  path = c(
    "ipw_data/easy_setting_pubcov_ipw_test",
    "ipw_data/medium_setting_pubcov_ipw_test",
    "ipw_data/hard_setting_pubcov_ipw_test"
  )
)

combine_intermediate_files <- function(path, pattern = "intermediate_result_.*\\.RData") {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No intermediate files found in: ", path)
  }

  bind_rows(lapply(files, function(f) {
    e <- new.env()
    load(f, envir = e)
    e[[ls(e)[1]]]
  }))
}

results_df <- bind_rows(lapply(seq_len(nrow(setting_dirs)), function(i) {
  combine_intermediate_files(setting_dirs$path[i]) %>%
    mutate(setting = setting_dirs$setting[i])
})) %>%
  mutate(
    setting = factor(setting, levels = c("Favorable", "Typical", "Extreme")),
    PUMA = as.character(PUMA)
  )

acs_pop <- read_csv("data/ACS_NPS_pop.csv", show_col_types = FALSE) %>%
  mutate(PUMA = as.character(PUMA))

true_values <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(
    response_true = mean(.data[[response_var]], na.rm = TRUE),
    .groups = "drop"
  )

metrics <- true_values %>%
  left_join(results_df, by = "PUMA") %>%
  group_by(setting, model) %>%
  summarize(
    n_sims = n_distinct(sim_num),
    Bias = mean(point_est - response_true, na.rm = TRUE),
    Abs_Bias = mean(abs(point_est - response_true), na.rm = TRUE),
    MSE = mean((response_true - point_est)^2, na.rm = TRUE),
    MAB = mean(abs(response_true - point_est), na.rm = TRUE),
    Coverage = mean(between(response_true, lower_CI, upper_CI), na.rm = TRUE),
    IS = mean(int_score(alpha, response_true, lower_CI, upper_CI), na.rm = TRUE),
    .groups = "drop"
  )

metrics <- metrics %>%
  mutate(
    family = case_when(
      model == "direst" ~ "Direct",
      str_detect(model, "E&V RE") ~ "E&V RE",
      str_detect(model, "E&V") ~ "E&V",
      str_detect(model, "ALP") ~ "ALP",
      str_detect(model, "CLW") ~ "CLW",
      str_detect(model, "calibration") ~ "Calibration",
      TRUE ~ "Other"
    ),
    label = case_when(
      model == "direst" ~ "Direct",

      model == "IPW (E&V)" ~ "IPW-DE",
      model == "IPW (E&V)+Y" ~ "IPW-DE+Y",
      model == "IPW (E&V RE)" ~ "IPW-DE (RE)",
      model == "IPW (E&V RE)+Y" ~ "IPW-DE+Y (RE)",

      model == "IPW (ALP) [AGEP]" ~ "ALP AGEP",
      model == "IPW (CLW) [AGEP]" ~ "CLW AGEP",
      model == "IPW (calibration) [AGEP]" ~ "CAL AGEP",

      model == "IPW (ALP) [AGEP+Y]" ~ "ALP AGEP+Y",
      model == "IPW (CLW) [AGEP+Y]" ~ "CLW AGEP+Y",
      model == "IPW (calibration) [AGEP+Y]" ~ "CAL AGEP+Y",

      model == "IPW (ALP) [AGEP_binned]" ~ "ALP AGEP-bin",
      model == "IPW (CLW) [AGEP_binned]" ~ "CLW AGEP-bin",
      model == "IPW (calibration) [AGEP_binned]" ~ "CAL AGEP-bin",

      model == "IPW (ALP) [AGEP_binned+Y]" ~ "ALP AGEP-bin+Y",
      model == "IPW (CLW) [AGEP_binned+Y]" ~ "CLW AGEP-bin+Y",
      model == "IPW (calibration) [AGEP_binned+Y]" ~ "CAL AGEP-bin+Y",

      TRUE ~ model
    ),
    family = factor(
      family,
      levels = c("Direct", "ALP", "CLW", "Calibration", "E&V", "E&V RE", "Other")
    )
  )

family_colors <- c(
  "Direct" = "grey35",
  "ALP" = "#0072B2",
  "CLW" = "#009E73",
  "Calibration" = "#CC79A7",
  "E&V" = "#E69F00",
  "E&V RE" = "#D55E00",
  "Other" = "black"
)

family_shapes <- c(
  "Direct" = 4,
  "ALP" = 16,
  "CLW" = 17,
  "Calibration" = 15,
  "E&V" = 18,
  "E&V RE" = 8,
  "Other" = 1
)

base_theme <- theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.text = element_text(color = "black"),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )

baselines <- metrics %>%
  filter(model == "direst") %>%
  mutate(ref_label = "Direct")

plot_df <- metrics %>%
  filter(model != "direst")

method_order <- plot_df %>%
  filter(setting == "Typical") %>%
  arrange(family, MSE) %>%
  pull(label) %>%
  unique()

if (length(method_order) == 0) {
  method_order <- plot_df %>%
    arrange(family, MSE) %>%
    pull(label) %>%
    unique()
}

add_direct_line <- function(p, bl_df, x_var) {
  bl_df <- bl_df %>% rename(.x = !!sym(x_var))

  p +
    geom_vline(
      data = bl_df,
      aes(xintercept = .x),
      color = "grey35",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    geom_text(
      data = bl_df,
      aes(x = .x, y = Inf, label = ref_label),
      color = "grey35",
      size = 3,
      angle = 90,
      vjust = -0.5,
      hjust = 1.25,
      inherit.aes = FALSE
    )
}

make_metric_plot <- function(metric, x_label, filename, log_scale = FALSE, add_zero = FALSE) {
  p <- plot_df %>%
    filter(!is.na(.data[[metric]]), !is.nan(.data[[metric]])) %>%
    mutate(label = factor(label, levels = rev(method_order))) %>%
    ggplot(aes(x = .data[[metric]], y = label, color = family, shape = family)) +
    geom_point(size = 3) +
    facet_wrap(~ setting, scales = "free_x") +
    scale_color_manual(values = family_colors, name = "Method") +
    scale_shape_manual(values = family_shapes, name = "Method") +
    labs(x = x_label, y = NULL) +
    coord_cartesian(clip = "off") +
    base_theme

  if (log_scale) {
    p <- p + scale_x_log10()
  }

  if (add_zero) {
    p <- p + geom_vline(xintercept = 0, color = "red", linewidth = 0.5)
  }

  bl_df <- baselines %>%
    filter(!is.na(.data[[metric]]), !is.nan(.data[[metric]]))

  p <- add_direct_line(p, bl_df, metric)

  ggsave(
    filename,
    p,
    width = 12,
    height = 7,
    dpi = 600
  )

  p
}

p_mse <- make_metric_plot(
  metric = "MSE",
  x_label = "MSE (log scale)",
  filename = "figures/fig_ipw_test_mse_dotplot.pdf",
  log_scale = TRUE
)

p_bias <- make_metric_plot(
  metric = "Bias",
  x_label = "Bias",
  filename = "figures/fig_ipw_test_bias_dotplot.pdf",
  add_zero = TRUE
)

p_abs_bias <- make_metric_plot(
  metric = "Abs_Bias",
  x_label = "Mean absolute bias",
  filename = "figures/fig_ipw_test_abs_bias_dotplot.pdf"
)

p_is <- make_metric_plot(
  metric = "IS",
  x_label = "Interval Score (log scale)",
  filename = "figures/fig_ipw_test_is_dotplot.pdf",
  log_scale = TRUE
)

p_cov <- plot_df %>%
  filter(!is.na(Coverage), !is.nan(Coverage)) %>%
  mutate(label = factor(label, levels = rev(method_order))) %>%
  ggplot(aes(x = Coverage, y = label, color = family, shape = family)) +
  geom_vline(
    data = baselines %>% filter(!is.na(Coverage), !is.nan(Coverage)),
    aes(xintercept = Coverage),
    color = "grey35",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  geom_text(
    data = baselines %>% filter(!is.na(Coverage), !is.nan(Coverage)),
    aes(x = Coverage, y = Inf, label = ref_label),
    color = "grey35",
    size = 3,
    angle = 90,
    vjust = -0.5,
    hjust = 1.25,
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = 0.95,
    color = "red",
    linewidth = 0.5
  ) +
  geom_point(size = 3) +
  facet_wrap(~ setting, scales = "free_x") +
  scale_color_manual(values = family_colors, name = "Method") +
  scale_shape_manual(values = family_shapes, name = "Method") +
  labs(x = "Coverage", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme

ggsave(
  "figures/fig_ipw_test_coverage_dotplot.pdf",
  p_cov,
  width = 12,
  height = 7,
  dpi = 600
)
