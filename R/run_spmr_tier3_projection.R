#!/usr/bin/env Rscript

# Generate all seven Tier 3 SPM alternatives from a saved RTMB pollock fit.
# The installed spmR package supplies the public run and table interfaces.

suppressPackageStartupMessages(library(spmR))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

rtmb_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(rtmb_root, "R", "Rpm.R"))) {
  stop("Run this script from the RTMB pollock repository root.")
}

writer_env <- new.env(parent = globalenv())
sys.source(file.path(rtmb_root, "R", "write_spmr_projection_inputs.R"),
           envir = writer_env)

model_file <- Sys.getenv(
  "POLLOCK_RTMB_RESULT",
  file.path(rtmb_root, "analysis", "output", "base.rds")
)
output_dir <- Sys.getenv(
  "POLLOCK_SPMR_OUTPUT",
  file.path(rtmb_root, "analysis", "output", "spmR_projection")
)

handoff <- writer_env$write_spmr_projection_inputs(
  model_file = model_file,
  output_dir = output_dir,
  config_path = file.path(rtmb_root, "R", "config.R"),
  alt_list = 1:7
)

detail <- spmR::runSPM(handoff$output_dir, run = TRUE, engine = "admb")
projection_years <- sort(unique(detail$Year))
first_comparison <- which(projection_years > max(handoff$fixed_catch_years))[1]
comparison_years <- projection_years[
  first_comparison:min(length(projection_years), first_comparison + 1L)
]

tier3_table <- spmR::tier3_scenario_table(
  detail,
  years = comparison_years,
  digits = 2
)
readr::write_csv(
  tier3_table,
  file.path(handoff$output_dir, "tier3_seven_scenario_table.csv")
)

print(tier3_table, n = Inf, width = Inf)
message("Tier 3 table written to ", file.path(
  handoff$output_dir, "tier3_seven_scenario_table.csv"
))

# Recreate the separately displayed Alternative 2 fixed-catch trajectory from
# the same refreshed bridge output.
alt2_output_dir <- file.path(
  rtmb_root, "analysis", "output", "spmR_projection_alt2_fixed1300"
)
alt2_fixed_catches <- stats::setNames(rep(1300, length(2025:2032)), 2025:2032)
alt2_handoff <- writer_env$write_spmr_projection_inputs(
  model_file = model_file,
  output_dir = alt2_output_dir,
  config_path = file.path(rtmb_root, "R", "config.R"),
  alt_list = 2L,
  fixed_catches = alt2_fixed_catches,
  nproj_years = length(alt2_fixed_catches),
  run_name = "rtmb_alt2_fixed1300"
)
alt2_detail <- spmR::runSPM(alt2_handoff$output_dir, run = TRUE, engine = "admb")
readr::write_csv(
  alt2_detail,
  file.path(alt2_handoff$output_dir, "spm_detail.csv")
)
message("Alternative 2 fixed-catch detail written to ", file.path(
  alt2_handoff$output_dir, "spm_detail.csv"
))
