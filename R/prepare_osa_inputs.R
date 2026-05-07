source(file.path("R", "config.R"))

model_file <- file.path(rtmb_dir, "analysis", "output", "base.rds")
output_file <- file.path(rtmb_dir, "analysis", "output", "osa_inputs.rds")

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

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(osa_inputs, output_file)
message("Wrote OSA inputs to ", output_file)
