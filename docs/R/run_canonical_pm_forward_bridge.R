#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
output_dir <- Sys.getenv(
  "RCEATTLE_BRIDGE_OUTPUT_DIR",
  file.path(project_root, "results", "canonical_pm", "bridge_forward_pass")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  RCEATTLE_LIB = file.path(project_root, ".r-lib-rceattle-5.8.1"),
  RCEATTLE_CANONICAL_PM = "false",
  RCEATTLE_ADMB_DIR = file.path(project_root, "ADMB", "m23_rceattle_full")
)

source(file.path(project_root, "2024 EBS pollock bridging.R"), chdir = FALSE)

likelihood_components <- rowSums(bridging_model_1$quantities$jnll_comp)
admb_report_lines <- readLines(file.path(AD, "pm.rep"), warn = FALSE)
admb_parameter_header <- readLines(
  file.path(AD, "pm.par"),
  n = 1L,
  warn = FALSE
)
read_admb_scalar <- function(key) {
  key_row <- which(admb_report_lines == key)[1]
  if (is.na(key_row) || key_row == length(admb_report_lines)) {
    return(NA_real_)
  }
  as.numeric(strsplit(trimws(admb_report_lines[key_row + 1L]), " +")[[1]][1])
}
dummy_gradient <- tryCatch(
  max(abs(bridging_model_1$obj$gr())),
  error = function(e) NA_real_
)
admb_maximum_gradient <- as.numeric(sub(
  ".*Maximum gradient component = ([^ ]+).*$",
  "\\1",
  admb_parameter_header
))

comparison <- data.frame(
  Metric = c(
    "N ratio minimum",
    "N ratio maximum",
    "SSB ratio minimum",
    "SSB ratio maximum",
    "SSB mean absolute percent difference",
    "Catch mean absolute percent difference",
    "Rceattle fixed-pass objective (all scientific parameters fixed)",
    "Rceattle reported index negative log-likelihood",
    "Rceattle reported catch negative log-likelihood",
    "Rceattle reported composition negative log-likelihood",
    "Rceattle reported recruitment-deviation negative log-likelihood",
    "Rceattle reported likelihood-component sum",
    "ADMB total negative log-likelihood (tot_like; not directly comparable)",
    "ADMB data negative log-likelihood (dat_like; not directly comparable)",
    "ADMB maximum gradient component",
    "Rceattle dummy-parameter maximum absolute gradient (not a convergence diagnostic)",
    "Estimated scientific parameters in fixed pass"
  ),
  Value = c(
    min(N_rce / N_admb),
    max(N_rce / N_admb),
    min(ssb_rce / ssb_admb),
    max(ssb_rce / ssb_admb),
    100 * mean(abs(ssb_rce / ssb_admb - 1)),
    100 * mean(abs(cat_rce / pred_cat - 1)),
    as.numeric(bridging_model_1$quantities$jnll),
    likelihood_components[["Index data"]],
    likelihood_components[["Catch data"]],
    likelihood_components[["Composition data"]],
    likelihood_components[["Recruitment deviates"]],
    sum(likelihood_components),
    read_admb_scalar("tot_like"),
    read_admb_scalar("dat_like"),
    admb_maximum_gradient,
    dummy_gradient,
    0
  )
)

trajectory_comparison <- data.frame(
  Year = yrs,
  SSB_ADMB = ssb_admb,
  SSB_Rceattle = ssb_rce,
  SSB_ratio = ssb_rce / ssb_admb,
  Catch_ADMB = pred_cat,
  Catch_Rceattle = cat_rce,
  Catch_ratio = cat_rce / pred_cat
)

saveRDS(
  list(
    model = bridging_model_1,
    data = fp,
    comparison = comparison,
    trajectories = trajectory_comparison,
    N_ADMB = N_admb,
    N_Rceattle = N_rce,
    selectivity = list(fishery = sel_fsh, BTS = sel_bts, ATS = sel_ats),
    admb_reference = AD
  ),
  file.path(output_dir, "bridge_forward_pass_results.rds")
)
write.csv(
  comparison,
  file.path(output_dir, "bridge_forward_pass_summary.csv"),
  row.names = FALSE
)
write.csv(
  trajectory_comparison,
  file.path(output_dir, "bridge_forward_pass_trajectories.csv"),
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "bridge_forward_pass_session_info.txt")
)
