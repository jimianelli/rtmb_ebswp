#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
rceattle_lib <- Sys.getenv(
  "RCEATTLE_LIB",
  file.path(project_root, ".r-lib-rceattle-5.8.1")
)
.libPaths(c(rceattle_lib, .libPaths()))
suppressPackageStartupMessages(library(Rceattle))
stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))

output_dir <- file.path("results", "rceattle_5.8.1_validation", "self_test")
fit <- readRDS(file.path(output_dir, "upstream_afsc_base_fit.rds"))$fitted
initial <- readRDS(file.path(output_dir, "self_test_50_all.rds"))

message("Repeating the 50 simulations from the estimated-parameter start.")
all_sims <- list()
all_converged <- logical()
for (batch in seq_len(5L)) {
  batch_seed <- 20240814L + (batch - 1L) * 10L
  message("Restart comparison batch ", batch, " of 5.")
  sims <- self_test(
    fit,
    nsim = 10,
    seed = batch_seed,
    cores = 4,
    getsd = FALSE,
    phase = TRUE,
    start = "estimated",
    debug = TRUE,
    timeout = 1800
  )
  converged <- attr(sims, "converged")
  names(sims) <- sprintf("Sim_%02d", seq_along(sims) + (batch - 1L) * 10L)
  names(converged) <- names(sims)
  saveRDS(
    list(models = sims, converged = converged, seed_base = batch_seed),
    file.path(output_dir, sprintf("self_test_estimated_start_batch_%02d.rds", batch))
  )
  all_sims <- c(all_sims, sims)
  all_converged <- c(all_converged, converged)
}

objective <- function(model) {
  if (inherits(model, "condition")) return(NA_real_)
  value <- model$opt$objective
  if (is.null(value)) value <- model$opt$opt$objective
  as.numeric(value)
}
initial_objective <- vapply(initial$models, objective, numeric(1))
estimated_objective <- vapply(all_sims, objective, numeric(1))
comparison <- data.frame(
  Simulation = names(all_sims),
  Initial_start_converged = unname(initial$converged),
  Estimated_start_converged = unname(all_converged),
  Initial_start_objective = initial_objective,
  Estimated_start_objective = estimated_objective,
  Estimated_minus_initial_objective = estimated_objective - initial_objective,
  stringsAsFactors = FALSE
)
write.csv(
  comparison,
  file.path(output_dir, "self_test_restart_comparison.csv"),
  row.names = FALSE
)
saveRDS(
  list(models = all_sims, converged = all_converged),
  file.path(output_dir, "self_test_50_estimated_start_all.rds")
)
message("Paired restart comparison saved under ", output_dir, ".")
