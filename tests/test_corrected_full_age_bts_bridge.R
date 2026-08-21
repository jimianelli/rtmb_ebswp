#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd(), mustWork = TRUE)
run_dir <- file.path(
  repo_root, "analysis", "output", "corrected_full_age_bts",
  "admb_root", "runs", "full_age_bts"
)
rtmb_file <- file.path(
  repo_root, "analysis", "output", "corrected_full_age_bts", "rtmb_base.rds"
)

if (!file.exists(file.path(run_dir, "pm.rep")) || !file.exists(rtmb_file)) {
  stop("Run the corrected ADMB and RTMB bridge workflows first.")
}

fit <- readRDS(rtmb_file)
metrics <- fit$metadata$bridge_metrics
stopifnot(
  identical(fit$metadata$configuration, "corrected_full_age_bts"),
  identical(fit$metadata$bts_comp_normalization, "full_ages"),
  metrics$absolute_total_nll_difference < 0.02,
  metrics$maximum_key_percent_difference < 0.01,
  metrics$maximum_absolute_gradient < 1e-3,
  metrics$maximum_bts_observed_sum_error < 1e-10,
  metrics$maximum_bts_predicted_sum_error < 1e-10
)

cat("PASS: corrected full-age BTS ADMB-to-RTMB bridge reproduced.\n")
