#!/usr/bin/env Rscript

# Run SparseNUTS for the accepted base RTMB model. SparseNUTS controls all
# sampler settings; the only additional argument is the RTMB global data object
# needed to rebuild the objective in parallel R sessions.

suppressPackageStartupMessages(library(SparseNUTS))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
  error = function(e) NA_character_
)
rtmb_root <- if (!is.na(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

rtmb_env <- new.env(parent = globalenv())
rtmb_env$rm <- function(...) invisible(NULL)
rtmb_env$source <- function(file, ...) {
  base::source(file, local = parent.frame(), ...)
}
source(file.path(rtmb_root, "R", "config.R"), local = rtmb_env)

obj <- rtmb_env$obj
model_data <- rtmb_env$data
if (!identical(as.integer(model_data$fishery_sel_form %||% 0L), 0L)) {
  stop("The SparseNUTS default runner must use fishery_sel_form = 0 (base RTMB model).")
}

output_file <- file.path(
  rtmb_root, "analysis", "output", "sparsenuts",
  "rtmb_ebswp_sparsenuts_default.rds"
)
monitor_file <- file.path(rtmb_root, "reporting", "data", "sparsenuts_base_monitor.csv")
summary_file <- file.path(rtmb_root, "reporting", "data", "sparsenuts_base_run_summary.csv")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(monitor_file), showWarnings = FALSE, recursive = TRUE)

package_description <- utils::packageDescription("SparseNUTS")
cat(
  "Running the base RTMB model with SparseNUTS package defaults.\n",
  "Package version: ", as.character(utils::packageVersion("SparseNUTS")), "\n",
  "Package commit: ", package_description$RemoteSha %||% "not recorded", "\n",
  "Parallel RTMB globals: model_data exported as 'data'.\n",
  sep = ""
)

started <- Sys.time()
fit <- SparseNUTS::sample_snuts(
  obj,
  globals = list(data = model_data),
  model_name = "EBS pollock base RTMB"
)
finished <- Sys.time()

attr(fit, "rtmb_ebswp_sparsenuts") <- list(
  runner = "R/run_sparsenuts_default.R",
  model = "accepted base RTMB model (fishery_sel_form = 0)",
  fishery_sel_form = 0L,
  package_version = as.character(utils::packageVersion("SparseNUTS")),
  package_remote_sha = package_description$RemoteSha %||% NA_character_,
  call = paste(
    "SparseNUTS::sample_snuts(obj,",
    "globals = list(data = model_data),",
    "model_name = 'EBS pollock base RTMB')"
  ),
  sampler_settings = "SparseNUTS package defaults",
  execution = "parallel chains using package-default chains and cores",
  globals = "data",
  started = started,
  finished = finished,
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs"))
)

tmp_file <- paste0(output_file, ".tmp")
saveRDS(fit, tmp_file)
if (!file.rename(tmp_file, output_file)) {
  stop("Could not move completed SparseNUTS fit into place: ", output_file)
}
diagnostics <- SparseNUTS::check_snuts_diagnostics(fit, print = FALSE)
utils::write.csv(as.data.frame(fit$monitor), monitor_file, row.names = FALSE)
utils::write.csv(
  data.frame(
    model = "accepted base RTMB model (fishery_sel_form = 0)",
    package_version = as.character(utils::packageVersion("SparseNUTS")),
    package_remote_sha = package_description$RemoteSha %||% NA_character_,
    algorithm = fit$algorithm,
    metric = fit$metric,
    chains = dim(fit$samples)[2],
    warmup_per_chain = fit$warmup,
    post_warmup_per_chain = dim(fit$samples)[1] - fit$warmup,
    divergences = round(
      dim(fit$samples)[2] * (dim(fit$samples)[1] - fit$warmup) *
        diagnostics$perc_divergent[1] / 100
    ),
    divergence_percent = diagnostics$perc_divergent[1],
    maximum_rhat = max(fit$monitor$rhat, na.rm = TRUE),
    minimum_bulk_ess = min(fit$monitor$ess_bulk, na.rm = TRUE),
    median_bulk_ess = median(fit$monitor$ess_bulk, na.rm = TRUE),
    minimum_tail_ess = min(fit$monitor$ess_tail, na.rm = TRUE),
    elapsed_seconds = attr(fit, "rtmb_ebswp_sparsenuts")$elapsed_seconds
  ),
  summary_file,
  row.names = FALSE
)
cat("Saved completed SparseNUTS fit: ", output_file, "\n", sep = "")
cat("Saved current monitor and run summary tables in reporting/data.\n")
