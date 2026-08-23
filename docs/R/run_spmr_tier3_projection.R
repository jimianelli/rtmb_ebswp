#!/usr/bin/env Rscript

.libPaths(c(
  file.path(getwd(), ".r-lib-rceattle-5.8.1"),
  file.path(getwd(), ".r-lib-pm"),
  .libPaths()
))
suppressPackageStartupMessages(library(spmR))
source(file.path("R", "write_spmr_projection_inputs.R"))

main <- write_spmr_projection_inputs()
detail <- spmR::runSPM(main$output_dir, run = TRUE, engine = "admb")
readr::write_csv(detail, file.path(main$output_dir, "spm_detail.csv"))
projection_means <- detail |>
  dplyr::group_by(Alt, Year) |>
  dplyr::summarize(
    B = mean(SSB),
    Catch = mean(Catch),
    ABC = mean(ABC),
    OFL = mean(OFL),
    F = mean(F),
    B_B35_percent = 100 * mean(SSB / B35),
    .groups = "drop"
  )
readr::write_csv(
  projection_means,
  file.path(main$output_dir, "spm_projection_means.csv")
)
years <- sort(unique(detail$Year))
comparison_years <- head(years[years > max(main$fixed_catch_years)], 2)
tier3 <- spmR::tier3_scenario_table(detail, years = comparison_years, digits = 2)
readr::write_csv(tier3, file.path(main$output_dir, "tier3_seven_scenario_table.csv"))

fixed <- stats::setNames(rep(1300, length(2025:2032)), 2025:2032)
alt2 <- write_spmr_projection_inputs(
  output_dir = "results/canonical_pm/spmR_projection_alt2_fixed1300", alt_list = 2L,
  fixed_catches = fixed, nproj_years = length(fixed), run_name = "rceattle_alt2_fixed1300"
)
alt2_detail <- spmR::runSPM(alt2$output_dir, run = TRUE, engine = "admb")
readr::write_csv(alt2_detail, file.path(alt2$output_dir, "spm_detail.csv"))

print(tier3, n = Inf, width = Inf)
