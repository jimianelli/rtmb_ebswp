#!/usr/bin/env Rscript

# Rebuild the corrected full-age BTS ADMB-to-RTMB bridge and every primary
# downstream product used by the RTMB working paper.

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "R", "config.R"))) {
  stop("Run this script from the rtmb_ebswp repository root.")
}

rscript <- file.path(R.home("bin"), "Rscript")
case_dir <- file.path(repo_root, "analysis", "output", "corrected_full_age_bts")
admb_run_dir <- file.path(case_dir, "admb_root", "runs", "full_age_bts")
base_file <- file.path(case_dir, "rtmb_base.rds")
only_bts_file <- file.path(case_dir, "only_bts.rds")
retro_file <- file.path(case_dir, "retro_9_peel.rds")
osa_input_file <- file.path(case_dir, "osa_inputs.rds")
osa_output_dir <- file.path(case_dir, "osa")
osa_output_file <- file.path(osa_output_dir, "rtmb_ebswp_osa_residuals.rds")
projection_dir <- file.path(case_dir, "spmR_projection")
projection_alt2_dir <- file.path(case_dir, "spmR_projection_alt2_fixed1300")

run_step <- function(label, script, env = character()) {
  message("\n== ", label, " ==")
  status <- system2(
    rscript,
    file.path(repo_root, script),
    env = env,
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop(label, " failed with status ", status)
  }
}

common_env <- c(
  "EBSWP_BRIDGE_CASE=corrected_full_age_bts",
  paste0("EBSWP_ADMB_RUN_DIR=", admb_run_dir),
  "EBSWP_BTS_COMP_NORMALIZATION=full_ages"
)

skip_base <- tolower(Sys.getenv("RTMB_SKIP_CORRECTED_BASE", unset = "false")) %in%
  c("1", "true", "yes")
if (!skip_base) {
  run_step(
    "Corrected full-age BTS ADMB model",
    "scripts/run_corrected_full_age_bts_admb.R"
  )
  run_step(
    "Corrected ADMB-to-RTMB bridge",
    "R/write_output.R",
    c(common_env, "EBSWP_OUTPUT_FILE=corrected_full_age_bts/rtmb_base.rds")
  )
}
run_step(
  "Corrected bridge regression gate",
  "tests/test_corrected_full_age_bts_bridge.R"
)
run_step(
  "Only-BTS sensitivity",
  "R/run_only_bts.R",
  c(
    common_env,
    paste0("RTMB_ONLY_BTS_OUTPUT=", only_bts_file),
    paste0("RTMB_ONLY_BTS_BASE_FILE=", base_file)
  )
)
run_step(
  "OSA input preparation",
  "R/prepare_osa_inputs.R",
  c(
    common_env,
    paste0("RTMB_OSA_MODEL_FILE=", base_file),
    paste0("RTMB_OSA_INPUT_FILE=", osa_input_file)
  )
)
run_step(
  "OSA residual diagnostics",
  "R/run_osa_comps.R",
  c(
    paste0("RTMB_OSA_INPUT_FILE=", osa_input_file),
    paste0("RTMB_OSA_OUTPUT_DIR=", osa_output_dir),
    paste0("RTMB_OSA_OUTPUT_FILE=", osa_output_file)
  )
)
run_step(
  "Nine-peel retrospective",
  "R/run_retrospective.R",
  c(
    common_env,
    paste0("RTMB_RETRO_BASE_FILE=", base_file),
    paste0("RTMB_RETRO_OUTPUT=", retro_file),
    "RTMB_RETRO_PEELS=0:9",
    "RTMB_RETRO_STOP_ON_ERROR=true"
  )
)
run_step(
  "Seven-scenario and fixed-catch spmR projections",
  "R/run_spmr_tier3_projection.R",
  c(
    common_env,
    paste0("POLLOCK_RTMB_RESULT=", base_file),
    paste0("POLLOCK_SPMR_OUTPUT=", projection_dir),
    paste0("POLLOCK_SPMR_ALT2_OUTPUT=", projection_alt2_dir)
  )
)

base_md5 <- unname(tools::md5sum(base_file))
retro <- readRDS(retro_file)
osa <- readRDS(osa_output_file)
projection_lineage <- read.csv(
  file.path(projection_dir, "base_lineage.csv"),
  stringsAsFactors = FALSE
)
alt2_lineage <- read.csv(
  file.path(projection_alt2_dir, "base_lineage.csv"),
  stringsAsFactors = FALSE
)
only_bts <- readRDS(only_bts_file)

lineage <- data.frame(
  product = c(
    "Corrected RTMB bridge", "Only-BTS sensitivity", "OSA diagnostics",
    "Nine-peel retrospective", "Seven-scenario projections",
    "Alternative 2 fixed-catch projections"
  ),
  base_md5 = c(
    base_md5,
    only_bts$metadata$base_md5,
    osa$lineage$base_md5,
    retro$base_lineage$md5,
    projection_lineage$model_md5[1],
    alt2_lineage$model_md5[1]
  ),
  stringsAsFactors = FALSE
)
lineage$matches_corrected_base <- lineage$base_md5 == base_md5
if (any(!lineage$matches_corrected_base)) {
  stop("A downstream product has a different corrected-base checksum.")
}
write.csv(lineage, file.path(case_dir, "downstream_lineage.csv"), row.names = FALSE)

if (tolower(Sys.getenv("RTMB_RENDER_REPORT", unset = "true")) %in%
    c("1", "true", "yes")) {
  message("\n== Render corrected RTMB working paper ==")
  status <- system2(
    "quarto",
    c("render", file.path("reporting", "ebs_pollock_rtmb_ebswp_assessment.qmd")),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop("Quarto render failed with status ", status)
  }
}

message("\nCorrected full-age BTS rebuild completed.")
message("Base MD5: ", base_md5)
message("Lineage: ", file.path(case_dir, "downstream_lineage.csv"))
