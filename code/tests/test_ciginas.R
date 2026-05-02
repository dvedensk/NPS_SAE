library(readr)
library(dplyr)
library(survey)
library(sampling)

source(file.path("code", "sampling_functions.R"))
source(file.path("code", "utils.R"))
source(file.path("code", "models", "ciginas.R"))

set.seed(99)

response_var <- "PUBCOV"
X_formula <- as.formula("~ AGEP_binned + RAC1P + SEX")

PS_weight_config <- list(WAGP = 0.05, PWGTP = -0.2)
NPS_weight_config <- list(PWGTP = 0.10, POVPIP = -1.52)

acs_pop <- read_csv(file.path("data", "ACS_NPS_pop.csv")) %>%
  mutate(
    AGEP_binned = factor(AGEP_binned),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

truth <- acs_pop %>%
  group_by(PUMA) %>%
  summarize(truth = mean(.data[[response_var]], na.rm = TRUE), .groups = "drop")

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

ps_df <- acs_pop[ps_sample$idx, ] %>%
  mutate(weights = ps_sample$weights)

nps_df <- acs_pop[nps_sample$idx, ]

ciginas_res <- ciginas_out(
  ps = ps_df,
  nps = nps_df,
  X_formula = X_formula,
  response = response_var,
  ps_weight_var = "weights"
)

required_cols <- c(
  "PUMA", "ps_est", "nps_ipw_est", "weight_nps", "point_est",
  "lower_CI", "upper_CI", "model"
)
stopifnot(all(required_cols %in% names(ciginas_res)))
stopifnot(nrow(ciginas_res) > 0)
stopifnot(!anyDuplicated(ciginas_res$PUMA))
stopifnot(all(is.finite(stats::na.omit(ciginas_res$point_est))))
stopifnot(all(stats::na.omit(ciginas_res$weight_nps) >= 0))
stopifnot(all(stats::na.omit(ciginas_res$weight_nps) <= 1))
stopifnot(all(stats::na.omit(ciginas_res$point_est) >= 0))
stopifnot(all(stats::na.omit(ciginas_res$point_est) <= 1))

convex_rows <- ciginas_res %>%
  filter(!is.na(ps_est), !is.na(nps_ipw_est))

stopifnot(all(
  convex_rows$point_est >= pmin(convex_rows$ps_est, convex_rows$nps_ipw_est) - 1e-8
))
stopifnot(all(
  convex_rows$point_est <= pmax(convex_rows$ps_est, convex_rows$nps_ipw_est) + 1e-8
))

metric_summary <- ciginas_res %>%
  left_join(truth, by = "PUMA") %>%
  summarize(
    ps_mae = mean(abs(ps_est - truth), na.rm = TRUE),
    nps_ipw_mae = mean(abs(nps_ipw_est - truth), na.rm = TRUE),
    ciginas_mae = mean(abs(point_est - truth), na.rm = TRUE)
  )

print(metric_summary)
print(ciginas_res %>% select(PUMA, ps_est, nps_ipw_est, weight_nps, point_est) %>% head())
cat("\nCiginas smoke test passed.\n")
