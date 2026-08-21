#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), ".r-lib-rceattle-5.8.1"), .libPaths()))
suppressPackageStartupMessages(library(Rceattle))
stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

output_dir <- file.path("results", "rceattle_5.8.1_validation", "extended_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fits <- readRDS(file.path("results", "canonical_pm", "ebs_pollock_method_fits.rds"))
fit <- fits$nonparametric_pm
cores <- max(1L, min(8L, parallel::detectCores() - 2L))

message("Running 50 phased jitters with Rceattle 5.8.1.")
set.seed(20260819)
jitters <- Rceattle::jitter(
  Rceattle = fit,
  njitter = 50,
  phase = TRUE,
  seed = 20260819,
  cores = cores,
  getsd = FALSE,
  timeout = 1800
)
saveRDS(jitters, file.path(output_dir, "jitter_50.rds"))

jitter_models <- jitters$Rceattle_list
if (is.null(jitter_models)) jitter_models <- list()
jitter_nll <- jitters$nll
if (is.null(jitter_nll)) jitter_nll <- numeric()
jitter_summary <- data.frame(
  Run = seq_along(jitter_nll),
  Objective = as.numeric(jitter_nll),
  Delta_objective = as.numeric(jitter_nll) - min(as.numeric(jitter_nll), na.rm = TRUE),
  Convergence_status = if (length(jitter_models)) {
    vapply(
      jitter_models,
      function(model) model$convergence$status %||% NA_character_,
      character(1)
    )
  } else {
    rep(NA_character_, length(jitter_nll))
  }
)
write.csv(jitter_summary, file.path(output_dir, "jitter_50_summary.csv"), row.names = FALSE)

message("Profiling age-3+ natural mortality with Rceattle 5.8.1.")
profile_m <- stats::profile(
  fitted = fit,
  param = "M1",
  slots = list(c(1, 1, 3)),
  values = list(seq(0.15, 0.50, by = 0.05)),
  cores = cores,
  getsd = FALSE
)
saveRDS(profile_m, file.path(output_dir, "profile_M_age3plus.rds"))
write.csv(
  transform(
    as.data.frame(profile_m$grid),
    Objective = as.numeric(profile_m$nll),
    Delta_objective = as.numeric(profile_m$nll) -
      min(as.numeric(profile_m$nll), na.rm = TRUE)
  ),
  file.path(output_dir, "profile_M_age3plus.csv"),
  row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
message("Extended diagnostics saved under ", output_dir, ".")
