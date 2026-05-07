# Write projection inputs for the spmR/SPM workflow from saved RTMB output.
# Usage: Rscript R/write_spmr_projection_inputs.R

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

format_spm_numbers <- function(x) {
  paste(format(signif(as.numeric(x), 8), scientific = FALSE, trim = TRUE), collapse = " ")
}

write_pm_prj <- function(report, data, output_file, stock_name = NULL) {
  stopifnot(!is.null(report$N), !is.null(report$F), !is.null(report$SSB))

  terminal_year <- as.integer(data$endyr)
  nages <- as.integer(data$nages)
  if (is.null(stock_name)) {
    stock_name <- paste0("EBS_Pollock_", terminal_year)
  }

  terminal_index <- nrow(report$N)
  terminal_rows <- max(1, terminal_index - 4):terminal_index
  f_age <- min(6, ncol(report$F))
  avg_f <- mean(as.numeric(report$F[terminal_rows, f_age]), na.rm = TRUE)

  wt_fut <- data$wt_fut %||% data$wt_fsh[terminal_index, ]
  sel_fut <- as.numeric(report$sel_fsh[terminal_index, ])
  natage <- as.numeric(report$N[terminal_index, ])
  rec <- as.numeric(report$N[, 1])
  ssb <- as.numeric(report$SSB)

  rec_years <- data$styr:data$endyr
  rec_keep <- rec_years >= 1978 & rec_years <= terminal_year
  ssb_keep <- rec_years >= 1977 & rec_years <= terminal_year - 1

  lines <- c(
    stock_name,
    "1    # SSLn species...",
    "0    # Buffer of Dorn",
    "1    # Number of fsheries",
    "1    # Number of sexes",
    paste0(format_spm_numbers(avg_f), "  # averagei 5yr f"),
    "1  # author f",
    "0.4  # ABC SPR",
    "0.35 # MSY/OFL SPR",
    "4  # Spawnmo",
    paste0(nages, " # Number of ages"),
    "1  # Fratio",
    paste0(format_spm_numbers(data$natmort), " # Natural Mortality"),
    "# Maturity",
    format_spm_numbers(data$p_mature / max(data$p_mature, na.rm = TRUE)),
    "# Wt spawn",
    format_spm_numbers(data$wt_ssb[terminal_index, ]),
    "# Wt fsh",
    format_spm_numbers(wt_fut),
    "# selectivity",
    format_spm_numbers(sel_fut),
    "# natage",
    format_spm_numbers(natage),
    "# Nrec",
    as.character(sum(rec_keep)),
    "# rec",
    format_spm_numbers(rec[rec_keep]),
    "# SpawningBiomass",
    format_spm_numbers(ssb[ssb_keep])
  )

  writeLines(lines, output_file)
  invisible(output_file)
}

write_spm_dat <- function(output_file,
                          spp_file = "pm.prj",
                          begin_year,
                          alt_list = 1:7,
                          fixed_catches = c(1350, 1350),
                          nproj_years = 14,
                          nsims = 1000,
                          run_name = "rtmb_ebswp") {
  if (is.null(names(fixed_catches))) {
    catch_years <- begin_year + seq_along(fixed_catches) - 1L
  } else {
    catch_years <- as.integer(names(fixed_catches))
  }
  fixed_catches <- as.numeric(fixed_catches)
  alt_list <- as.integer(alt_list)

  lines <- c(
    "#_SETUP_FILE_FOR_EBS_Pollock_RTMB",
    paste(run_name, "# Run name"),
    "# Tier",
    "3",
    "#--------------------------------------------",
    paste(length(alt_list), "# Number of Alternatives"),
    "#--------------------------------------------",
    paste(alt_list, "# List of alternatives"),
    "#--------------------------------------------",
    "1    # Flag to set TAC equal to ABC",
    "#--------------------------------------------",
    "2    # Stock-recruitment type (1=Ricker, 2=Bholt)",
    "1    # projection recruitment form",
    "1    # SR-Conditioning",
    "0    # Recruitment prior CV condition",
    "#--------------------------------------------",
    "1    # Flag to write big file",
    "#--------------------------------------------",
    paste0(nproj_years, " #_Number of projection years"),
    paste0(nsims, " #_Number of simulations"),
    paste0(begin_year, " #_Begin Year"),
    "#_Number_of_years with specified catch",
    as.character(length(fixed_catches)),
    "# Number of species",
    "1",
    "# OY Min",
    "0",
    "# OY Max",
    "2.00E+06",
    "# data files for each species",
    spp_file,
    "# ABC Multipliers",
    "1",
    "# scalars",
    "1",
    "# New Alt 4 Fabc SPRs",
    "0.6",
    "# Number of TAC model categories",
    "1",
    "# TAC model indices",
    "1",
    "# Catch in each future year",
    paste(catch_years, fixed_catches)
  )

  writeLines(lines, output_file)
  invisible(output_file)
}

write_spmr_projection_inputs <- function(model_file = file.path("analysis", "output", "base.rds"),
                                         output_dir = file.path("analysis", "output", "spmR_projection"),
                                         config_path = NULL,
                                         template_dir = NULL,
                                         alt_list = 1:7,
                                         fixed_catches = c(1350, 1350),
                                         nproj_years = 14,
                                         nsims = 1000,
                                         run_name = "rtmb_ebswp") {
  if (is.null(config_path)) {
    config_path <- file.path("R", "config.R")
  }
  config_path <- normalizePath(config_path, mustWork = TRUE)

  env <- new.env(parent = globalenv())
  env$rm <- function(...) invisible(NULL)
  env$source <- function(file, ...) {
    base::source(file, local = parent.frame(), ...)
  }
  source(config_path, local = env)

  rtmb_dir <- env$rtmb_dir
  data <- env$data

  model_path <- normalizePath(file.path(rtmb_dir, model_file), mustWork = TRUE)
  output_dir <- normalizePath(file.path(rtmb_dir, output_dir), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  saved <- readRDS(model_path)
  report <- saved$report %||% saved$rtmb %||% saved

  pm_prj <- file.path(output_dir, "pm.prj")
  spm_dat <- file.path(output_dir, "spm.dat")
  tacpar <- file.path(output_dir, "tacpar.dat")
  spm_exe <- file.path(output_dir, "spm")

  write_pm_prj(report = report, data = data, output_file = pm_prj)
  write_spm_dat(
    output_file = spm_dat,
    begin_year = as.integer(data$endyr) + 1L,
    alt_list = alt_list,
    fixed_catches = fixed_catches,
    nproj_years = nproj_years,
    nsims = nsims,
    run_name = run_name
  )

  if (is.null(template_dir)) {
    template_dir <- file.path(env$pollock_root, "admb", "runs", "for_rtmb", "proj")
  }
  if (file.exists(file.path(template_dir, "tacpar.dat"))) {
    file.copy(file.path(template_dir, "tacpar.dat"), tacpar, overwrite = TRUE)
  }
  if (file.exists(file.path(template_dir, "spm"))) {
    file.copy(file.path(template_dir, "spm"), spm_exe, overwrite = TRUE)
    Sys.chmod(spm_exe, mode = "0755")
  }

  manifest <- data.frame(
    file = c(spm_dat, pm_prj, tacpar, spm_exe),
    role = c(
      "SPM setup file",
      "Species assessment input file generated from RTMB",
      "SPM TAC parameter file copied from ADMB bridge projection directory",
      "SPM executable copied from ADMB bridge projection directory"
    ),
    exists = file.exists(c(spm_dat, pm_prj, tacpar, spm_exe)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE)

  invisible(list(output_dir = output_dir, manifest = manifest))
}

is_rscript_call <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  any(grepl("^--file=.*write_spmr_projection_inputs\\.R$", args))
}

if (!interactive() && is_rscript_call()) {
  result <- write_spmr_projection_inputs()
  message("Wrote spmR projection inputs to ", result$output_dir)
}
