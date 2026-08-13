library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(knitr)
library(kableExtra)

source(file.path("code", "utils.R"))

alpha <- 0.05
response_var <- "PUBCOV"

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

combine_intermediate_files <- function(path, pattern = "intermediate_result_.*\\.RData", setting) {
  files <- list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No intermediate result files found in: ", path)
  }

  data_list <- lapply(files, function(f) {
    e <- new.env()
    load(f, envir = e)
    out <- e[[ls(e)[1]]]

    if (!"sim_num" %in% names(out)) {
      sim_num <- as.integer(str_match(basename(f), "intermediate_result_([0-9]+)\\.RData")[, 2])
      out$sim_num <- sim_num
    }

    out
  })

  bind_rows(data_list) %>%
    mutate(setting = setting)
}

easy_df <- combine_intermediate_files(
  path = "results/easy_setting_pubcov",
  setting = "Favorable"
)

medium_df <- combine_intermediate_files(
  path = "results/medium_setting_pubcov",
  setting = "Typical"
)

hard_df <- combine_intermediate_files(
  path = "results/hard_setting_pubcov",
  setting = "Extreme"
)

results_df <- bind_rows(easy_df, medium_df, hard_df) %>%
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

model_map <- tibble(
  model = c(
    "Direct",
    "BULM-PS",
    "IPW-BULM",
    "IPW-BULM+Y",
    "IPW-DE",
    "IPW-DE+Y",
    "IPW-DE-CLIP",
    "MRP-INT-P",
    "MRP-INT-P-IPW",
    "NIP (exp link)",
    "NIP (p-value)",
    "VSW"
  ),
  #same as model names now, but can edit labels below if necessary
  #while upstream names above will be unchanged
  label = c(
    "Direct",
    "BULM-PS",
    "IPW-BULM",
    "IPW-BULM+Y",
    "IPW-DE",
    "IPW-DE+Y",
    "IPW-DE-CLIP",
    "MRP-INT-P",
    "MRP-INT-P-IPW",
    "NIP (exp link)",
    "NIP (p-value)",
    "VSW"
  ),
  family = c(
    "Direct",
    "BULM",
    "BULM",
    "BULM",
    "IPW",
    "IPW",
    "IPW",
    "MRP",
    "MRP",
    "NIP",
    "NIP",
    "VSW"
  ),
  role = c(
    "reference",
    "reference",
    "method",
    "method",
    "method",
    "method",
    "method",
    "method",
    "method",
    "method",
    "method",
    "method"
  )
)

metrics <- metrics %>%
  inner_join(model_map, by = "model") %>%
  mutate(
    family = factor(
      family,
      levels = c("Direct", "BULM", "IPW", "MRP", "NIP", "VSW")
    ),
    setting = factor(
      setting,
      levels = c("Favorable", "Typical", "Extreme")
    )
  )

family_order <- c("BULM", "IPW", "MRP", "NIP", "VSW")

family_colors <- c(
  "Direct" = "grey40",
  "BULM" = "#E69F00",
  "IPW" = "#56B4E9",
  "MRP" = "#009E73",
  "NIP" = "#CC79A7",
  "VSW" = "#D55E00"
)

family_shapes <- c(
  "Direct" = 4,
  "BULM" = 16,
  "IPW" = 17,
  "MRP" = 15,
  "NIP" = 18,
  "VSW" = 8
)

base_theme <- theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
  )

baselines <- metrics %>%
  filter(role == "reference") %>%
  mutate(
    ref_label = recode(
      label,
      "Direct" = "D.Est",
      "BULM-PS" = "BULM-PS"
    )
  )

plot_df <- metrics %>%
  filter(role == "method")

method_order <- plot_df %>%
  filter(setting == "Typical", family %in% family_order) %>%
  mutate(family = factor(family, levels = family_order)) %>%
  arrange(family, MSE) %>%
  pull(label) %>%
  unique()

if (length(method_order) == 0) {
  method_order <- plot_df %>%
    mutate(family = factor(family, levels = family_order)) %>%
    arrange(family, MSE) %>%
    pull(label) %>%
    unique()
}

baseline_col <- "grey40"

add_baselines <- function(p, bl_df, x_var) {
  bl_df <- bl_df %>% rename(.x = !!sym(x_var))

  p +
    geom_vline(
      data = bl_df,
      aes(xintercept = .x),
      color = baseline_col,
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    geom_text(
      data = bl_df,
      aes(x = .x, y = Inf, label = ref_label),
      color = baseline_col,
      size = 2.8,
      angle = 90,
      vjust = -0.5,
      hjust = 1.3,
      inherit.aes = FALSE
    )
}

p_mse <- plot_df %>%
  mutate(label = factor(label, levels = rev(method_order))) %>%
  ggplot(aes(x = MSE, y = label, color = family, shape = family)) +
  geom_point(size = 2.5) +
  facet_wrap(~ setting, scales = "free_x") +
  scale_color_manual(values = family_colors[family_order], name = "Method Family") +
  scale_shape_manual(values = family_shapes[family_order], name = "Method Family") +
  scale_x_log10() +
  labs(x = "MSE (log scale)", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme

p_mse <- add_baselines(p_mse, baselines, "MSE")

ggsave(
  "figures/fig_mse_dotplot.pdf",
  p_mse,
  width = 12,
  height = 6,
  dpi = 600
)

baselines_is <- baselines %>%
  filter(!is.nan(IS), !is.na(IS))

p_is <- plot_df %>%
  filter(!is.nan(IS), !is.na(IS)) %>%
  mutate(label = factor(label, levels = rev(method_order))) %>%
  ggplot(aes(x = IS, y = label, color = family, shape = family)) +
  geom_point(size = 2.5) +
  facet_wrap(~ setting, scales = "free_x") +
  scale_color_manual(values = family_colors[family_order], name = "Method Family") +
  scale_shape_manual(values = family_shapes[family_order], name = "Method Family") +
  scale_x_log10() +
  labs(x = "Interval Score (log scale)", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme

p_is <- add_baselines(p_is, baselines_is, "IS")

ggsave(
  "figures/fig_is_dotplot.pdf",
  p_is,
  width = 12,
  height = 6,
  dpi = 600
)

baselines_cov <- baselines %>%
  filter(!is.nan(Coverage), !is.na(Coverage))

p_cov <- plot_df %>%
  filter(!is.nan(Coverage), !is.na(Coverage)) %>%
  mutate(label = factor(label, levels = rev(method_order))) %>%
  ggplot(aes(x = Coverage, y = label, color = family, shape = family)) +
  geom_vline(
    data = baselines_cov,
    aes(xintercept = Coverage),
    color = baseline_col,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  geom_text(
    data = baselines_cov,
    aes(x = Coverage, y = Inf, label = ref_label),
    color = baseline_col,
    size = 2.8,
    angle = 90,
    vjust = -0.5,
    hjust = 1.3,
    inherit.aes = FALSE
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ setting, scales = "free_x") +
  scale_color_manual(values = family_colors[family_order], name = "Method Family") +
  scale_shape_manual(values = family_shapes[family_order], name = "Method Family") +
  geom_vline(
    xintercept = 0.95,
    color = "red",
    linewidth = 0.5
  ) +
  labs(x = "Coverage", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme

ggsave(
  "figures/fig_coverage_dotplot.pdf",
  p_cov,
  width = 12,
  height = 6,
  dpi = 600
)

baselines_mab <- baselines %>%
  filter(!is.nan(MAB), !is.na(MAB))

p_mab <- plot_df %>%
  filter(!is.nan(MAB), !is.na(MAB)) %>%
  mutate(label = factor(label, levels = rev(method_order))) %>%
  ggplot(aes(x = MAB, y = label, color = family, shape = family)) +
  geom_vline(
    data = baselines_mab,
    aes(xintercept = MAB),
    color = baseline_col,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  geom_text(
    data = baselines_mab,
    aes(x = MAB, y = Inf, label = ref_label),
    color = baseline_col,
    size = 2.8,
    angle = 90,
    vjust = -0.5,
    hjust = 1.3,
    inherit.aes = FALSE
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ setting, scales = "free_x") +
  scale_color_manual(values = family_colors[family_order], name = "Method Family") +
  scale_shape_manual(values = family_shapes[family_order], name = "Method Family") +
  labs(x = "MAB", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme

ggsave(
  "figures/fig_bias_dotplot.pdf",
  p_mab,
  width = 12,
  height = 6,
  dpi = 600
)

write_csv(
  metrics %>% arrange(setting, MSE),
  "figures/metrics_main_clean.csv"
)

base_data <- metrics %>%
  mutate(
    Method = factor(
      label,
      levels = c(
        "Direct",
        "BULM-PS",
        "IPW-DE",
        "IPW-DE+Y",
        "IPW-DE-CLIP",
        "IPW-BULM",
        "IPW-BULM+Y",
        "MRP-INT-P-IPW",
        "MRP-INT-P",
        "NIP (exp link)",
        "NIP (p-value)",
        "VSW"
      )
    )
  ) %>%
  filter(!is.na(Method)) %>%
  arrange(Method)

n_rep <- max(base_data$n_sims, na.rm = TRUE)

tab1_formatted <- base_data %>%
  mutate(MSE = MSE * 1000) %>%
  select(Method, setting, MSE) %>%
  pivot_wider(names_from = setting, values_from = MSE) %>%
  mutate(
    Favorable_min = min(Favorable, na.rm = TRUE),
    Typical_min = min(Typical, na.rm = TRUE),
    Extreme_min = min(Extreme, na.rm = TRUE)
  ) %>%
  mutate(
    Favorable = sprintf("%.1f", Favorable),
    Typical = sprintf("%.1f", Typical),
    Extreme = sprintf("%.1f", Extreme),
    Favorable_min = sprintf("%.1f", Favorable_min),
    Typical_min = sprintf("%.1f", Typical_min),
    Extreme_min = sprintf("%.1f", Extreme_min)
  ) %>%
  mutate(
    Favorable = ifelse(Favorable == Favorable_min, paste0("\\textbf{", Favorable, "}"), Favorable),
    Typical = ifelse(Typical == Typical_min, paste0("\\textbf{", Typical, "}"), Typical),
    Extreme = ifelse(Extreme == Extreme_min, paste0("\\textbf{", Extreme, "}"), Extreme)
  ) %>%
  mutate(` ` = "", Method = as.character(Method)) %>%
  select(` `, Method, Favorable, Typical, Extreme)

tab1_latex <- tab1_formatted %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "llccc",
    col.names = c("", "Method", "Favorable", "Typical", "Extreme"),
    caption = paste0(
      "\\label{tab:mse}\n",
      "MSE ($\\times 10^3$) by DDC setting.\n",
      "Results from ", n_rep, " replications.\n",
      "Bold indicates lowest MSE within each setting.\n"
    ),
    linesep = ""
  ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 9) %>%
  add_header_above(c(" " = 2, "DDC Setting" = 3), line = FALSE) %>%
  row_spec(0, extra_latex_after = "\\cmidrule(l{3pt}r{3pt}){3-5}") %>%
  row_spec(c(2, 7, 9, 11), extra_latex_after = "\\midrule")

tab2_formatted <- base_data %>%
  select(Method, setting, Coverage, IS) %>%
  mutate(
    Cov_fmt = ifelse(is.na(Coverage), NA_character_, sub("^0", "", sprintf("%.3f", Coverage))),
    IS_fmt = ifelse(is.na(IS), NA_character_, sprintf("%.3f", IS))
  ) %>%
  group_by(setting) %>%
  mutate(min_is = min(IS, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    IS_final = case_when(
      is.na(IS) ~ NA_character_,
      round(IS, 3) == round(min_is, 3) ~ paste0("\\textbf{", sub("^0", "", IS_fmt), "}"),
      TRUE ~ sub("^0", "", IS_fmt)
    ),
    cell_val = case_when(
      is.na(Coverage) | is.na(IS) ~ "---",
      TRUE ~ paste0(Cov_fmt, " (", IS_final, ")")
    )
  ) %>%
  select(Method, setting, cell_val) %>%
  pivot_wider(names_from = setting, values_from = cell_val) %>%
  mutate(` ` = "", Method = as.character(Method)) %>%
  select(` `, Method, Favorable, Typical, Extreme)

tab2_latex <- tab2_formatted %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "llccc",
    col.names = c("", "Method", "Favorable", "Typical", "Extreme"),
    caption = paste0(
      "\\label{tab:interval}\n",
      "Coverage (Interval Score) by DDC setting.\n",
      "Results from ", n_rep, " replications.\n",
      "Nominal coverage is 95\\%.\n",
      "Bold indicates lowest IS within each setting."
    ),
    linesep = ""
  ) %>%
  kable_styling(latex_options = c("hold_position"), font_size = 9) %>%
  add_header_above(c(" " = 2, "DDC Setting" = 3), line = FALSE) %>%
  row_spec(0, extra_latex_after = "\\cmidrule(l{3pt}r{3pt}){3-5}") %>%
  row_spec(c(2, 7, 9, 11), extra_latex_after = "\\midrule")

cat(tab1_latex, file = "tables/table_mse.tex")
cat(tab2_latex, file = "tables/table_coverage.tex")

metrics %>%
  count(setting, model, label, family, role) %>%
  arrange(setting, family, label) %>%
  print(n = Inf)
