# ----------- LOADING DATA -----------
# Load required libraries and functions
library(readr)
library(dplyr)
library(cmdstanr)

# Load test sample data
load("data/test_sample.RData")

# Load population data for MRP
acs_pop <- read_csv("data/ACS_NPS_pop.csv") %>%
  mutate(
    AGEP = factor(AGEP),
    RAC1P = factor(RAC1P),
    SEX = factor(SEX),
    PUMA = factor(PUMA)
  )

# ----------- DESCRIPTIVE STATS -----------
# Show sample coverage by PUMA
cat("\nSample coverage:\n")
cat("PS covers", length(unique(ps_df$PUMA)), "out of", length(unique(acs_pop$PUMA)), "PUMAs\n")
cat("NPS covers", length(unique(nps_df$PUMA)), "out of", length(unique(acs_pop$PUMA)), "PUMAs\n")

# Show health insurance coverage rates
cat("\nHealth insurance coverage rates:\n")
cat("Population HICOV rate:", round(mean(acs_pop$HICOV), 3), "\n")
cat("PS HICOV rate:", round(mean(ps_df$HICOV), 3), "\n")
cat("NPS HICOV rate:", round(mean(nps_df$HICOV), 3), "\n")

cat("\n✅ Test data loaded and ready for MRP analysis!\n")

# ----------- MRP SET UP -----------
# Source MRP functions
source(file.path("code", "models", "mrp_all.R"))

# Load compiled Stan models
mod <- cmdstan_model(
  file.path("code", "models", "si2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

modINT <- cmdstan_model(
  file.path("code", "models", "mrp_int2.stan"),
  cpp_options = list(stan_threads = TRUE)
)

# Set up parallel processing
Sys.setenv(STAN_NUM_THREADS = parallel::detectCores())

# ----------- PROFILE MRP  -----------
cat("\n=== PROFILING MRP FUNCTIONS ===\n")

# Fix data objects: getMRP expects data frames, not list objects
nps_data <- nps_df # Use the data frame
ps_data <- ps_df # Use the data frame

cat("Data preparation:\n")
cat("- NPS data:", nrow(nps_data), "observations\n")
cat("- PS data:", nrow(ps_data), "observations\n")
cat("- Population data:", nrow(acs_pop), "observations\n\n")

# Profile getMRP function
cat("⏱️  Profiling getMRP()...\n")
start_time <- Sys.time()

tryCatch(
  {
    mrp <- getMRP(
      MR = nps_data,
      ps = ps_data,
      acs_pop = acs_pop,
      mod = mod
    )

    mrp_time <- Sys.time() - start_time
    cat("✅ getMRP completed in", round(mrp_time, 2), attr(mrp_time, "units"), "\n")

    # Post Stratification results (MRP-R & MRP-P)
    mrpr <- mrp$puma_summary_mrpr
    mrpp <- mrp$puma_summary_mrpp

    cat("- MRP-R results:", nrow(mrpr), "PUMAs\n")
    cat("- MRP-P results:", nrow(mrpp), "PUMAs\n")
  },
  error = function(e) {
    cat("❌ getMRP failed:", e$message, "\n")
  }
)

cat("\n⏱️  Profiling getMRP_INT()...\n")
start_time <- Sys.time()

tryCatch(
  {
    mrp1 <- getMRP_INT(
      MR = nps_data,
      ps = ps_data,
      acs_pop = acs_pop,
      mod = modINT
    )

    mrp_int_time <- Sys.time() - start_time
    cat("✅ getMRP_INT completed in", round(mrp_int_time, 2), attr(mrp_int_time, "units"), "\n")

    # Results
    mrpint <- mrp1$puma_summary_mrpp
    cat("- MRP-INT results:", nrow(mrpint), "PUMAs\n")

    # Show Rhat diagnostics if available
    if (!is.null(mrp1$rhat)) {
      max_rhat <- max(mrp1$rhat, na.rm = TRUE)
      cat("- Max R-hat:", round(max_rhat, 3), "\n")
      if (max_rhat > 1.1) {
        cat("⚠️  Warning: Some chains may not have converged (R-hat > 1.1)\n")
      }
    }
  },
  error = function(e) {
    cat("❌ getMRP_INT failed:", e$message, "\n")
  }
)

cat("\n⏱️  Profiling getMRP_INT(include_response = TRUE)...\n")
start_time <- Sys.time()

tryCatch(
  {
    mrp2 <- getMRP_INT(
      MR = nps_data,
      ps = ps_data,
      acs_pop = acs_pop,
      mod = modINT,
      include_response = TRUE
    )

    mrp_int_resp_time <- Sys.time() - start_time
    cat("✅ getMRP_INT(include_response=TRUE) completed in",
        round(mrp_int_resp_time, 2), attr(mrp_int_resp_time, "units"), "\n")

    mrpint_p <- mrp2$puma_summary_mrpp
    mrpint_r <- mrp2$puma_summary_mrpr
    cat("- MRP-INT-P (PUBCOV) results:", nrow(mrpint_p), "PUMAs\n")
    cat("- MRP-INT-R (PUBCOV) results:", nrow(mrpint_r), "PUMAs\n")
    cat("- Model labels:", unique(mrpint_p$model), "/", unique(mrpint_r$model), "\n")

    if (!is.null(mrp2$rhat)) {
      max_rhat <- max(mrp2$rhat, na.rm = TRUE)
      cat("- Max R-hat:", round(max_rhat, 3), "\n")
      if (max_rhat > 1.1) cat("⚠️  Warning: Some chains may not have converged (R-hat > 1.1)\n")
    }
  },
  error = function(e) {
    cat("❌ getMRP_INT(include_response=TRUE) failed:", e$message, "\n")
  }
)

cat("\n=== PROFILING COMPLETE ===\n")
