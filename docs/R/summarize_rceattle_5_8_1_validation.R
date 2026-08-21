#!/usr/bin/env Rscript

output_dir <- file.path("results", "rceattle_5.8.1_validation", "self_test")
base <- readRDS(file.path(output_dir, "upstream_afsc_base_fit.rds"))$fitted
initial <- readRDS(file.path(output_dir, "self_test_50_all.rds"))
restart <- read.csv(file.path(output_dir, "self_test_restart_comparison.csv"))
truncation <- read.csv(
  file.path(output_dir, "index_truncation_diagnostics_by_fleet.csv")
)

base_ssb <- as.numeric(base$quantities$ssb[1, ])
ssb_bias <- vapply(
  initial$models,
  function(model) {
    100 * median(as.numeric(model$quantities$ssb[1, ]) / base_ssb - 1)
  },
  numeric(1)
)
terminal_ssb_bias <- vapply(
  initial$models,
  function(model) {
    100 * (tail(as.numeric(model$quantities$ssb[1, ]), 1) /
      tail(base_ssb, 1) - 1)
  },
  numeric(1)
)

summary <- data.frame(
  Metric = c(
    "Initial-start self-tests converged",
    "Estimated-start self-tests converged",
    "Paired fits with objective difference greater than 0.0001",
    "Maximum absolute paired objective difference",
    "Median across-simulation median SSB percent bias",
    "Median terminal-year SSB percent bias",
    "Fleets with at least one row above 2 percent rejection probability",
    "Maximum row-level non-positive probability"
  ),
  Value = c(
    sum(initial$converged),
    sum(restart$Estimated_start_converged),
    sum(abs(restart$Estimated_minus_initial_objective) > 1e-4),
    max(abs(restart$Estimated_minus_initial_objective)),
    median(ssb_bias),
    median(terminal_ssb_bias),
    sum(truncation$Rows_above_2_percent > 0),
    max(truncation$Maximum_nonpositive_probability)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary,
  file.path(output_dir, "validation_summary.csv"),
  row.names = FALSE
)
