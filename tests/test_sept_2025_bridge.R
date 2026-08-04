#!/usr/bin/env Rscript

# Regression gate for the ADMB-to-RTMB bridge published in September 2025.
# Run from anywhere with:
#   Rscript tests/test_sept_2025_bridge.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this test with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
setwd(repo_root)

source(file.path("R", "config.R"))
data$return_nll_only <- FALSE
rtmb <- rpm(parms)

comparison <- compare_max_pct(rtmb, pm, tolerance = 1e-5)
comparison <- comparison[match(names(rtmb), comparison$variable), , drop = FALSE]

total_abs_difference <- abs(as.numeric(rtmb$tot_like) - as.numeric(pm$tot_like))
max_gradient <- max(abs(obj$gr()), na.rm = TRUE)

key_variables <- c(
  "N", "Z", "F", "M", "S", "pred_catch", "SSB", "pred_cpue",
  "pred_avo", "eb_bts", "eb_ats", "sel_fsh", "sel_bts", "sel_ats",
  "phat_fsh", "phat_bts", "phat_ats", "age_like", "cat_like",
  "bts_like", "ats_like", "rec_like", "sel_like", "sel_like_dev",
  "Priors", "tot_like"
)
key_comparison <- comparison[comparison$variable %in% key_variables, , drop = FALSE]
max_key_percent_difference <- max(key_comparison$max_abs_pct_diff, na.rm = TRUE)

cat(sprintf("ADMB total NLL: %.9f\n", as.numeric(pm$tot_like)))
cat(sprintf("RTMB total NLL: %.9f\n", as.numeric(rtmb$tot_like)))
cat(sprintf("Absolute total NLL difference: %.9f\n", total_abs_difference))
cat(sprintf("Maximum key percent difference: %.9f%%\n", max_key_percent_difference))
cat(sprintf("Maximum absolute gradient: %.9g\n", max_gradient))

failures <- character()
if (!is.finite(total_abs_difference) || total_abs_difference > 0.005) {
  failures <- c(failures, "total NLL difference exceeds 0.005")
}
if (!is.finite(max_key_percent_difference) || max_key_percent_difference > 0.0006) {
  failures <- c(failures, "a key output exceeds 0.0006 percent difference")
}
if (!is.finite(max_gradient) || max_gradient > 1e-5) {
  failures <- c(failures, "maximum absolute gradient exceeds 1e-5")
}
if (any(!key_comparison$equal)) {
  failures <- c(
    failures,
    paste0(
      "all.equal failed for: ",
      paste(key_comparison$variable[!key_comparison$equal], collapse = ", ")
    )
  )
}

output_path <- Sys.getenv("RTMB_BRIDGE_COMPARISON_CSV", unset = "")
if (nzchar(output_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(comparison, output_path, row.names = FALSE)
  cat("Wrote comparison:", normalizePath(output_path, mustWork = FALSE), "\n")
}

if (length(failures) > 0L) {
  stop("September 2025 bridge regression:\n- ", paste(failures, collapse = "\n- "))
}

cat("PASS: September 2025 ADMB-to-RTMB bridge reproduced.\n")
