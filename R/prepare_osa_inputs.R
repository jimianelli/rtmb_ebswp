args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
rtmb_dir <- if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  normalizePath(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

setwd(rtmb_dir)
source(file.path(rtmb_dir, "R", "config.R"))
rtmb_dir <- normalizePath(here::here(), mustWork = TRUE)

default_model_file <- if (identical(bridge_case, "corrected_full_age_bts")) {
  file.path(rtmb_dir, "analysis", "output", "corrected_full_age_bts", "rtmb_base.rds")
} else {
  file.path(rtmb_dir, "analysis", "output", "base.rds")
}
model_file <- Sys.getenv("RTMB_OSA_MODEL_FILE", unset = default_model_file)
output_file <- Sys.getenv(
  "RTMB_OSA_INPUT_FILE",
  unset = file.path(rtmb_dir, "analysis", "output", "osa_inputs.rds")
)

if (!file.exists(model_file)) {
  stop("Missing saved RTMB output at: ", model_file)
}

saved <- readRDS(model_file)
rtmb_report <- saved$report

osa_inputs <- list(
  Fishery = list(
    obs = data$oac_fsh,
    exp = rtmb_report$phat_fsh,
    n_eff = as.numeric(data$sam_fsh),
    index = seq_len(data$nages),
    years = data$yrs_fsh_data
  ),
  BTS = list(
    obs = data$oac_bts,
    exp = rtmb_report$phat_bts,
    n_eff = as.numeric(data$sam_bts),
    index = seq_len(data$nages),
    years = data$yrs_bts_data
  ),
  ATS = list(
    obs = data$oac_ats[, data$mina_ats:data$nages, drop = FALSE],
    exp = rtmb_report$phat_ats[, data$mina_ats:data$nages, drop = FALSE],
    n_eff = as.numeric(data$sam_ats),
    index = data$mina_ats:data$nages,
    years = data$yrs_ats_data
  )
)

# Exclude placeholder composition rows with no observed fish. afscOSA validates
# the rounded multinomial counts and correctly rejects these zero-information rows.
osa_keep_rows <- function(x) {
  totals <- rowSums(x$obs, na.rm = TRUE)
  rounded_counts <- round(sweep(x$obs, 1, x$n_eff / totals, `*`), 0)
  is.finite(x$n_eff) & x$n_eff >= 1 & totals > 0 & rowSums(rounded_counts) >= 1
}
dropped_rows <- lapply(osa_inputs, function(x) {
  keep <- osa_keep_rows(x)
  x$years[!keep]
})
osa_inputs <- lapply(osa_inputs, function(x) {
  keep <- osa_keep_rows(x)
  lapply(x, function(value) {
    if (is.matrix(value)) value[keep, , drop = FALSE] else if (length(value) == length(keep)) value[keep] else value
  })
})

attr(osa_inputs, "lineage") <- list(
  base_file = normalizePath(model_file, mustWork = TRUE),
  base_md5 = unname(tools::md5sum(model_file)),
  configuration = saved$metadata$configuration,
  bts_comp_normalization = saved$metadata$bts_comp_normalization,
  prepared = Sys.time(),
  r_version = R.version.string,
  excluded_zero_information_years = dropped_rows
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(osa_inputs, output_file)
message("Wrote OSA inputs to ", output_file)
