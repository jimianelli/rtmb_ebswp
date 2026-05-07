# Re-run the default SparseNUTS sampler for the RTMB-ADMB model.

suppressPackageStartupMessages({
  library(SparseNUTS)
})

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

output_file <- file.path(rtmb_root, "analysis", "output", "sparsenuts", "rtmb_ebswp_sparsenuts_default.rds")
figure_dir <- file.path(rtmb_root, "analysis", "output", "figures")
pairs_file <- file.path(figure_dir, "rtmb_ebswp_sparsenuts_pairs_slow.png")
marginals_file <- file.path(figure_dir, "rtmb_ebswp_sparsenuts_marginals_slow.png")

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

cat("Running SparseNUTS::sample_snuts(obj) with package defaults except cores = 1 for RTMB serial execution...\n")
fit <- SparseNUTS::sample_snuts(
  obj,
  cores = 1,
  globals = list(data = rtmb_env$data)
)
attr(fit, "rtmb_ebswp_sparsenuts") <- list(
  runner = "R/run_sparsenuts_default.R",
  r_version = R.version.string,
  call = "SparseNUTS::sample_snuts(obj, cores = 1, globals = list(data = data))",
  execution = "serial chains because RTMB parallel worker globals failed",
  created = Sys.time()
)

tmp_file <- paste0(output_file, ".tmp")
saveRDS(fit, tmp_file)
file.rename(tmp_file, output_file)
cat("Saved SparseNUTS fit:", output_file, "\n")

slow_names <- character()
slow_indices <- integer()
if (!is.null(fit$monitor) && !is.null(fit$samples)) {
  slow_names <- fit$monitor |>
    dplyr::arrange(dplyr::desc(rhat), ess_bulk) |>
    dplyr::filter(is.finite(rhat)) |>
    dplyr::slice_head(n = 12) |>
    dplyr::pull(variable)
  slow_indices <- match(slow_names, dimnames(fit$samples)[[3]])
  slow_indices <- slow_indices[is.finite(slow_indices)]
}

if (length(slow_indices) > 0) {
  pairs_indices <- slow_indices[seq_len(min(length(slow_indices), 6))]
  png(pairs_file, width = 1800, height = 1800, res = 180)
  pairs(fit, pars = pairs_indices, order = "slow", diag = "hist", plot = TRUE)
  dev.off()
  cat("Saved SparseNUTS pairs.tmbfit plot:", pairs_file, "\n")

  png(marginals_file, width = 1800, height = 1400, res = 180)
  SparseNUTS::plot_marginals(fit, pars = slow_indices, order = "slow", mfrow = c(3, 4))
  dev.off()
  cat("Saved SparseNUTS marginal plot:", marginals_file, "\n")
} else {
  warning("No SparseNUTS monitor diagnostics were available for slow-order plots.")
}
