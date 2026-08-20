#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
rceattle_lib <- Sys.getenv(
  "RCEATTLE_LIB",
  file.path(project_root, ".r-lib-rceattle-5.8.1")
)
.libPaths(c(rceattle_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Rceattle)
  library(dplyr)
})

stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))

output_dir <- file.path("results", "rceattle_5.8.1_validation", "self_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the final EBS fit produced with Rceattle 5.8.1.")
fits <- readRDS(file.path("results", "canonical_pm", "ebs_pollock_method_fits.rds"))
fit <- fits$nonparametric_pm
base <- list(stage_a = fits$stage_a, fitted = fit)
saveRDS(base, file.path(output_dir, "upstream_afsc_base_fit.rds"))

objective <- fit$opt$objective
if (is.null(objective)) objective <- fit$opt$opt$objective
max_gradient <- fit$convergence$checks$max_gradient$data$max_gradient
base_summary <- data.frame(
  Rceattle_version = as.character(packageVersion("Rceattle")),
  Objective = as.numeric(objective),
  Maximum_gradient = as.numeric(max_gradient),
  Convergence_status = fit$convergence$status,
  stringsAsFactors = FALSE
)
write.csv(
  base_summary,
  file.path(output_dir, "upstream_afsc_base_fit_summary.csv"),
  row.names = FALSE
)

# Quantify how often the fitted natural-scale normal marginals would generate
# invalid non-positive observations. This remains visible even though
# self_test() suppresses warnings raised by sim_mod() inside each worker.
idx <- fit$data_list$index_data
fc <- fit$data_list$fleet_control
idx$Index_distribution <- fc$Index_distribution[
  match(idx$Fleet_code, fc$Fleet_code)
]
idx$Prediction <- as.numeric(fit$quantities$index_hat)
idx$Reported_sd <- as.numeric(fit$quantities$log_index_sd)
idx$Fitted_row <- idx$Year > 0 & idx$Year <= fit$data_list$endyr &
  idx$Observation > 0
idx$Nonpositive_probability <- 0
normal_rows <- idx$Fitted_row & idx$Index_distribution == "Normal"
idx$Nonpositive_probability[normal_rows] <- pnorm(
  0,
  mean = idx$Prediction[normal_rows],
  sd = idx$Reported_sd[normal_rows]
)
mvn_fleets <- unique(idx$Fleet_code[
  idx$Fitted_row & idx$Index_distribution %in% c("MVN", "MVNORM")
])
for (fleet_code in mvn_fleets) {
  rows <- which(idx$Fitted_row & idx$Fleet_code == fleet_code)
  fleet_name <- idx$Fleet_name[rows[1]]
  sigma <- fit$data_list$index_cov[[fleet_name]]
  if (is.null(sigma)) {
    sigma <- fit$data_list$index_cov[[as.character(fleet_code)]]
  }
  if (!is.null(sigma) && nrow(sigma) == length(rows)) {
    idx$Reported_sd[rows] <- sqrt(diag(sigma))
    idx$Nonpositive_probability[rows] <- pnorm(
      0,
      mean = idx$Prediction[rows],
      sd = idx$Reported_sd[rows]
    )
  }
}

truncation_by_fleet <- idx |>
  filter(Fitted_row) |>
  group_by(Fleet_name, Fleet_code, Index_distribution) |>
  summarise(
    Rows = n(),
    Maximum_nonpositive_probability = max(Nonpositive_probability),
    Mean_nonpositive_probability = mean(Nonpositive_probability),
    Rows_above_2_percent = sum(Nonpositive_probability > 0.02),
    .groups = "drop"
  )
write.csv(
  idx,
  file.path(output_dir, "index_truncation_diagnostics_by_row.csv"),
  row.names = FALSE
)
write.csv(
  truncation_by_fleet,
  file.path(output_dir, "index_truncation_diagnostics_by_fleet.csv"),
  row.names = FALSE
)

message("Running 50 phased self-tests in five checkpointed batches.")
all_sims <- list()
all_converged <- logical()
for (batch in seq_len(5L)) {
  batch_seed <- 20240814L + (batch - 1L) * 10L
  message("Self-test batch ", batch, " of 5; seed base ", batch_seed, ".")
  sims <- self_test(
    fit,
    nsim = 10,
    seed = batch_seed,
    cores = 4,
    getsd = FALSE,
    phase = TRUE,
    start = "initial",
    debug = TRUE,
    timeout = 1800
  )
  converged <- attr(sims, "converged")
  names(sims) <- sprintf("Sim_%02d", seq_along(sims) + (batch - 1L) * 10L)
  names(converged) <- names(sims)
  saveRDS(
    list(models = sims, converged = converged, seed_base = batch_seed),
    file.path(output_dir, sprintf("self_test_batch_%02d.rds", batch))
  )
  all_sims <- c(all_sims, sims)
  all_converged <- c(all_converged, converged)
}

saveRDS(
  list(models = all_sims, converged = all_converged),
  file.path(output_dir, "self_test_50_all.rds")
)
self_test_summary <- data.frame(
  Simulation = names(all_sims),
  Converged = unname(all_converged),
  Error = vapply(
    all_sims,
    function(x) if (inherits(x, "condition")) conditionMessage(x) else NA_character_,
    character(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  self_test_summary,
  file.path(output_dir, "self_test_50_summary.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "session_info.txt")
)
message("Self-test validation saved under ", output_dir, ".")
