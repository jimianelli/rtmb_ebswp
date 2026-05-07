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

output_path <- file.path(output_dir, "base.rds")
saveRDS(
  list(
    report = rtmb_report,
    metadata = list(
      model = "rtmb_ebswp",
      created = Sys.time(),
      admb_rep = admb_rep_path,
      admb_par = admb_par_path
    )
  ),
  output_path
)

cat("Wrote RTMB-ADMB output:", output_path, "\n")
