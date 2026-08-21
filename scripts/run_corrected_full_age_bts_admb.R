#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "R", "config.R"))) {
  stop("Run this script from the rtmb_ebswp repository root.")
}

legacy_run <- file.path(repo_root, "admb", "runs", "for_rtmb")
legacy_data <- file.path(repo_root, "admb", "runs", "data")
output_root <- file.path(
  repo_root, "analysis", "output", "corrected_full_age_bts", "admb_root"
)
run_dir <- file.path(output_root, "runs", "full_age_bts")
data_dir <- file.path(output_root, "runs", "data")

if (dir.exists(output_root)) {
  unlink(output_root, recursive = TRUE)
}
dir.create(run_dir, recursive = TRUE)
dir.create(data_dir, recursive = TRUE)

run_inputs <- c("pm.tpl", "pm.dat", "pm.par", "control.dat", "compweights.ctl")
copied_run <- file.copy(
  file.path(legacy_run, run_inputs),
  file.path(run_dir, run_inputs),
  overwrite = TRUE
)
data_inputs <- list.files(legacy_data, full.names = TRUE)
copied_data <- file.copy(
  data_inputs,
  file.path(data_dir, basename(data_inputs)),
  overwrite = TRUE
)
if (!all(copied_run) || !all(copied_data)) {
  stop("Failed to stage the corrected ADMB run inputs.")
}

template_path <- file.path(run_dir, "pm.tpl")
legacy_template_md5 <- unname(tools::md5sum(template_path))
template <- readLines(template_path)
old_line <- "        oac_bts(i )     = oac_bts_data(i)/ot_bts(i);"
line_index <- which(template == old_line)
if (length(line_index) != 1L) {
  stop("The legacy BTS normalization line changed; review pm.tpl manually.")
}
replacement <- c(
  "        // Full-age BTS composition correction: observations and expectations",
  "        // both sum to one across ages 1--15. ot_bts remains the ages 2--15",
  "        // abundance-total observation used elsewhere in the model.",
  "        oac_bts(i )     = oac_bts_data(i)/sum(oac_bts_data(i));"
)
template <- append(template[-line_index], replacement, after = line_index - 1L)
writeLines(template, template_path)

run_command <- function(command, args = character()) {
  status <- system2(command, args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop(command, " failed with status ", status)
  }
}

old_wd <- setwd(run_dir)
on.exit(setwd(old_wd), add = TRUE)
run_command("admb", "pm.tpl")
run_command("./pm", c("-nox", "-iprint", "150", "-ainp", "pm.par"))
setwd(old_wd)

required_outputs <- file.path(run_dir, c("pm.par", "pm.rep", "pm.std"))
if (!all(file.exists(required_outputs))) {
  stop("The corrected ADMB run failed to create its required outputs.")
}

manifest <- data.frame(
  case = "corrected_full_age_bts",
  normalization = "observed and predicted BTS ages 1--15 sum to one",
  legacy_template = file.path(legacy_run, "pm.tpl"),
  legacy_template_md5 = legacy_template_md5,
  corrected_template = template_path,
  corrected_template_md5 = unname(tools::md5sum(template_path)),
  corrected_par_md5 = unname(tools::md5sum(file.path(run_dir, "pm.par"))),
  corrected_rep_md5 = unname(tools::md5sum(file.path(run_dir, "pm.rep"))),
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(dirname(output_root), "admb_corrected_manifest.csv"),
  row.names = FALSE
)

message("Corrected ADMB run: ", run_dir)
message("Corrected report MD5: ", manifest$corrected_rep_md5)
