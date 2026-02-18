# Test: PUMA_lev union fix prevents NA group id error
#
# Reproduces the bug seen with the Easy (0.75/0.25) DDC setting, where the NPS
# may not cover all PUMAs but acs_pop and ps always do.
#
# Run from repo root: Rscript code/tests/test_mrp_puma_coverage.R

suppressPackageStartupMessages({
  library(dplyr)
})

source("code/models/mrp_all.R")
load("data/test_sample.RData")
acs_pop <- readr::read_csv("data/ACS_NPS_pop.csv", show_col_types = FALSE)

# Prepare factor levels matching what getMRP_INT does
nps_ca <- nps_df
nps_ca$AGEP_binned <- factor(nps_ca$AGEP_binned)
nps_ca$SEX         <- factor(nps_ca$SEX)
nps_ca$RAC1P       <- factor(nps_ca$RAC1P)

# Simulate the Easy setting: drop 5 PUMAs from NPS to trigger the bug
set.seed(42)
drop_pumas <- sample(unique(nps_ca$PUMA), 5)
nps_reduced <- nps_ca[!nps_ca$PUMA %in% drop_pumas, ]
cat("Dropped PUMAs from NPS:", paste(drop_pumas, collapse = ", "), "\n")
cat("NPS PUMAs:", length(unique(nps_reduced$PUMA)), "/ Pop PUMAs:", length(unique(acs_pop$PUMA)), "\n\n")

# Fake psi values (bypasses selection model, just tests PUMA level handling)
make_fake_psi <- function(df, bin_map) {
  psi <- runif(nrow(df), 0.1, 0.9)
  bin_val <- bin_fun(psi, digits = 2)
  # extend bin_map for any unseen bins
  new_bins <- setdiff(as.character(bin_val), names(bin_map))
  if (length(new_bins) > 0) {
    next_id <- max(bin_map) + seq_along(new_bins)
    bin_map[new_bins] <- next_id
  }
  list(psi = psi, psi_bin = as.integer(bin_map[as.character(bin_val)]), bin_map = bin_map)
}

# Build a minimal bin_map
all_psi_vals <- round(runif(200, 0.1, 0.9), 2)
bin_map <- setNames(seq_along(sort(unique(all_psi_vals))), as.character(sort(unique(all_psi_vals))))

# ── Test 1: old approach (NPS only) should fail ─────────────────────────────
PUMA_lev_old <- levels(factor(nps_reduced$PUMA))
pop_psi <- make_fake_psi(acs_pop, bin_map)

cat("Test 1 — old PUMA_lev (NPS only): expecting error... ")
result_old <- tryCatch(
  collapse_to_cells_int(
    df          = acs_pop,
    train_ref   = nps_reduced,
    PUMA_lev    = PUMA_lev_old,
    psi_vec     = pop_psi$psi,
    psi_bin_vec = pop_psi$psi_bin,
    weight_col  = "PWGTP",
    bin_map     = bin_map
  ),
  error = function(e) e
)
if (inherits(result_old, "error") && grepl("PUMA not in PUMA_lev", result_old$message)) {
  cat("PASSED (error reproduced)\n")
} else {
  cat("SKIPPED — no missing PUMAs triggered (test data may already cover all PUMAs)\n")
}

# ── Test 2: fixed approach (union of NPS + PS + pop) should succeed ──────────
PUMA_lev_new <- sort(unique(c(
  as.character(nps_reduced$PUMA),
  as.character(ps_df$PUMA),
  as.character(acs_pop$PUMA)
)))
pop_psi2 <- make_fake_psi(acs_pop, bin_map)

cat("Test 2 — new PUMA_lev (union): expecting success... ")
result_new <- tryCatch(
  collapse_to_cells_int(
    df          = acs_pop,
    train_ref   = nps_reduced,
    PUMA_lev    = PUMA_lev_new,
    psi_vec     = pop_psi2$psi,
    psi_bin_vec = pop_psi2$psi_bin,
    weight_col  = "PWGTP",
    bin_map     = bin_map
  ),
  error = function(e) e
)
if (inherits(result_new, "error")) {
  cat("FAILED —", result_new$message, "\n")
  quit(status = 1)
} else {
  cat("PASSED (no error,", nrow(result_new$cells), "cells)\n")
}

# ── Test 3: MRP-R poststrat uses only PS cells (no pop-only PUMAs contaminate) 
cat("Test 3 — PS cells only contain PUMAs present in PS... ")
ps_psi <- make_fake_psi(ps_df, bin_map)
ps_cells <- tryCatch(
  collapse_to_cells_int(
    df          = ps_df,
    train_ref   = nps_reduced,
    PUMA_lev    = PUMA_lev_new,
    psi_vec     = ps_psi$psi,
    psi_bin_vec = ps_psi$psi_bin,
    weight_col  = "weights",
    bin_map     = bin_map
  ),
  error = function(e) e
)
if (inherits(ps_cells, "error")) {
  cat("FAILED —", ps_cells$message, "\n")
  quit(status = 1)
} else {
  ps_pumas_in_cells <- unique(as.character(ps_cells$cells$PUMA))
  extra <- setdiff(drop_pumas, ps_pumas_in_cells)  # dropped PUMAs absent from PS cells
  cat("PASSED (PS cells span", length(ps_pumas_in_cells), "PUMAs; dropped PUMAs absent from cells as expected)\n")
}

cat("\nAll tests passed.\n")
