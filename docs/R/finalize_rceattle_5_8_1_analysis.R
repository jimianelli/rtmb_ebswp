#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), ".r-lib-rceattle-5.8.1"), .libPaths()))

suppressPackageStartupMessages({
  library(Rceattle)
  library(dplyr)
})

stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))
source(file.path("R", "ebs_pollock_fit_helpers.R"))

output_dir <- file.path("results", "canonical_pm")
method_file <- file.path(output_dir, "ebs_pollock_method_fits.rds")
workbook <- file.path("Data", "EBS_24_pollock_m23_rceattle_full_1964-2024.xlsx")
admb_dir <- file.path("ADMB", "m23_rceattle_full")

fits <- readRDS(method_file)
base <- fits$nonparametric_pm

message("Fitting the standard-multinomial sensitivity with Rceattle 5.8.1.")
multinomial <- fit_ebs_pollock_base(
  xlsx = workbook,
  comp_loglike = "Multinomial",
  comp_offset = 0,
  verbose = 1
)

fits$data <- base$data_list
fits$standard_multinomial_stage_a <- multinomial$stage_a
fits$standard_multinomial <- multinomial$fitted
fits$admb_directory <- admb_dir
fits$years <- base$data_list$styr:base$data_list$endyr
fits$rceattle_version <- as.character(packageVersion("Rceattle"))
saveRDS(fits, method_file)

objective <- function(x) {
  value <- x$opt$objective
  if (is.null(value)) value <- x$opt$opt$objective
  if (!length(value)) return(NA_real_)
  as.numeric(value)
}

max_gradient <- function(x) {
  value <- x$convergence$checks$max_gradient$data$max_gradient
  if (!length(value)) return(NA_real_)
  as.numeric(value)
}

model_summary <- data.frame(
  Model = c(
    "Rceattle 5.8.1 final configuration",
    "Rceattle 5.8.1 standard multinomial sensitivity",
    "Rceattle 5.8.1 2D AR1 selectivity sensitivity"
  ),
  Composition_likelihood = c(
    "MultinomialAFSC",
    "Multinomial",
    "MultinomialAFSC"
  ),
  Offset = c(1e-3, 0, 1e-3),
  Objective = c(
    objective(base),
    objective(multinomial$fitted),
    objective(fits$ar1_2d)
  ),
  Maximum_gradient = c(
    max_gradient(base),
    max_gradient(multinomial$fitted),
    max_gradient(fits$ar1_2d)
  ),
  Convergence_status = c(
    base$convergence$status,
    multinomial$fitted$convergence$status,
    fits$ar1_2d$convergence$status
  )
)
write.csv(model_summary, file.path(output_dir, "model_summary.csv"), row.names = FALSE)

trajectory <- function(fit, label) {
  years <- fit$data_list$styr:fit$data_list$endyr
  data.frame(
    Model = label,
    Year = years,
    SSB = as.numeric(fit$quantities$ssb[1, seq_along(years)]),
    Recruitment = as.numeric(fit$quantities$R[1, seq_along(years)]),
    Biomass = as.numeric(fit$quantities$biomass[1, seq_along(years)])
  )
}

trajectories <- bind_rows(
  trajectory(base, "Rceattle 5.8.1 final configuration"),
  trajectory(multinomial$fitted, "Rceattle 5.8.1 standard multinomial"),
  trajectory(fits$ar1_2d, "Rceattle 5.8.1 2D AR1")
)
write.csv(trajectories, file.path(output_dir, "model_trajectories.csv"), row.names = FALSE)

package_source <- normalizePath("/Users/jim/_mymods/ceattle/Rceattle", mustWork = TRUE)
model_source <- normalizePath("/Users/jim/_mymods/ceattle/Rceattle-models", mustWork = TRUE)
rceattle_commit <- system2("git", c("-C", package_source, "rev-parse", "HEAD"), stdout = TRUE)
model_commit <- system2("git", c("-C", model_source, "rev-parse", "origin/master"), stdout = TRUE)
admb_report <- file.path(admb_dir, "pm.rep")

lineage <- data.frame(
  Item = c(
    "Rceattle version",
    "Rceattle upstream commit",
    "Rceattle source worktree",
    "Rceattle-models upstream commit",
    "Input workbook",
    "Input workbook MD5",
    "Modified ADMB report",
    "Modified ADMB report MD5"
  ),
  Value = c(
    as.character(packageVersion("Rceattle")),
    rceattle_commit,
    package_source,
    model_commit,
    normalizePath(workbook),
    unname(tools::md5sum(workbook)),
    normalizePath(admb_report),
    unname(tools::md5sum(admb_report))
  )
)
write.csv(lineage, file.path(output_dir, "lineage.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

message("Rceattle 5.8.1 analysis products finalized under ", output_dir, ".")
