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
