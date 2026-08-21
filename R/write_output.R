# Write standardized RTMB-ADMB output for bridging comparisons.
# Usage: Rscript R/write_output.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  config_path <- file.path(script_dir, "config.R")
} else {
  config_path <- file.path("R", "config.R")
}

source(config_path)

# config.R intentionally clears its sourcing environment to reproduce the
# historical bridge setup, so resolve the repository root again afterward.
rtmb_dir <- normalizePath(getwd(), mustWork = TRUE)
admb_rep_path <- normalizePath(file.path(admb_run_dir, "pm.rep"), mustWork = TRUE)
admb_par_path <- normalizePath(file.path(admb_run_dir, "pm.par"), mustWork = TRUE)

return_nll_only <- FALSE
data$return_nll_only <- 0
rtmb_result <- rpm(parms)

rtmb_report <- NULL
if (is.list(rtmb_result) && !is.null(rtmb_result$rtmb)) {
  rtmb_report <- rtmb_result$rtmb
} else if (is.list(rtmb_result) && !is.null(rtmb_result$report)) {
  rtmb_report <- rtmb_result$report
} else if (is.list(rtmb_result) && !is.null(rtmb_result$rep)) {
  rtmb_report <- rtmb_result$rep
} else if (is.list(rtmb_result) &&
           (!is.null(rtmb_result$SSB) || !is.null(rtmb_result$tot_like) || !is.null(rtmb_result$nll))) {
  rtmb_report <- rtmb_result
}

if (is.null(rtmb_report)) {
  detail <- if (is.list(rtmb_result)) paste(names(rtmb_result), collapse = ", ") else class(rtmb_result)
  stop("RTMB report is NULL. Ensure return_nll_only is FALSE and the model ran successfully. Returned: ", detail)
}

output_dir <- file.path(rtmb_dir, "analysis", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_name <- Sys.getenv(
  "EBSWP_OUTPUT_FILE",
  unset = if (bridge_case == "legacy_pm_bridge") {
    "base.rds"
  } else {
    file.path(bridge_case, "rtmb_base.rds")
  }
)
output_path <- file.path(output_dir, output_name)
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
bridge_comparison <- compare_max_pct(rtmb_report, pm, tolerance = 1e-5)
key_variables <- c(
  "N", "Z", "F", "M", "S", "pred_catch", "SSB", "pred_cpue",
  "pred_avo", "eb_bts", "eb_ats", "sel_fsh", "sel_bts", "sel_ats",
  "phat_fsh", "phat_bts", "phat_ats", "age_like", "cat_like",
  "bts_like", "ats_like", "rec_like", "sel_like", "sel_like_dev",
  "Priors", "tot_like"
)
key_comparison <- bridge_comparison[
  bridge_comparison$variable %in% key_variables, , drop = FALSE
]
bridge_metrics <- list(
  admb_total_nll = as.numeric(pm$tot_like),
  rtmb_total_nll = as.numeric(rtmb_report$tot_like),
  absolute_total_nll_difference = abs(
    as.numeric(rtmb_report$tot_like) - as.numeric(pm$tot_like)
  ),
  maximum_key_percent_difference = max(
    key_comparison$max_abs_pct_diff, na.rm = TRUE
  ),
  maximum_absolute_gradient = max(abs(obj$gr()), na.rm = TRUE),
  maximum_bts_observed_sum_error = max(
    abs(rowSums(data$oac_bts) - 1), na.rm = TRUE
  ),
  maximum_bts_predicted_sum_error = max(
    abs(rowSums(rtmb_report$phat_bts) - 1), na.rm = TRUE
  )
)
saveRDS(
  list(
    report = rtmb_report,
    bridge_comparison = bridge_comparison,
    metadata = list(
      model = "rtmb_ebswp",
      configuration = bridge_case,
      bts_comp_normalization = bts_comp_normalization,
      created = Sys.time(),
      admb_rep = admb_rep_path,
      admb_par = admb_par_path,
      bridge_metrics = bridge_metrics
    )
  ),
  output_path
)

cat("Wrote RTMB-ADMB output:", output_path, "\n")
cat(sprintf(
  "Bridge: |delta NLL| = %.9f; max key difference = %.9f%%; max gradient = %.9g\n",
  bridge_metrics$absolute_total_nll_difference,
  bridge_metrics$maximum_key_percent_difference,
  bridge_metrics$maximum_absolute_gradient
))
