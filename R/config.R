# =============================================================================
# EBS Pollock RTMB Configuration
# =============================================================================
# Configuration script for running the RTMB implementation
# Handles path setup and data loading for comparison with ADMB
# =============================================================================

rm(list = ls())

# --- Path Configuration ---
get_rtmb_root <- function() {
  env_root <- Sys.getenv("RTMB_EBSWP_ROOT", unset = NA_character_)
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  if (exists("config_path", inherits = TRUE)) {
    config_dir <- dirname(normalizePath(get("config_path", inherits = TRUE), mustWork = TRUE))
    return(normalizePath(file.path(config_dir, ".."), mustWork = TRUE))
  }

  script_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )
  if (!is.na(script_file)) {
    return(normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE))
  }

  wd <- normalizePath(getwd(), mustWork = TRUE)
  candidates <- c(wd, file.path(wd, "rtmb_ebswp"), file.path(wd, ".."))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "config.R"))) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }

  stop("Cannot locate this repository root. Set RTMB_EBSWP_ROOT.")
}

get_pollock_root <- function(rtmb_root) {
  env_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
  if (is.na(env_root) || !nzchar(env_root)) {
    env_root <- Sys.getenv("POLLOCK_BASE", unset = NA_character_)
  }
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  candidates <- c(
    dirname(rtmb_root),
    file.path(dirname(rtmb_root), "pollock"),
    file.path(rtmb_root, "data", "private", "pollock")
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "admb", "runs", "for_rtmb"))) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }

  stop(
    "Cannot locate the external pollock workspace. ",
    "Set POLLOCK_ROOT to a directory containing admb/runs/for_rtmb/."
  )
}

# Set up paths. This repository is standalone; the ADMB/input workspace is an
# external dependency supplied through POLLOCK_ROOT.
rtmb_dir <- get_rtmb_root()
pollock_root <- get_pollock_root(rtmb_dir)

cat("Project root:", pollock_root, "\n")
cat("RTMB directory:", rtmb_dir, "\n")

# --- Load Dependencies ---
library(RTMB)
library(tidyverse)
library(patchwork)

if (requireNamespace("TMBhelper", quietly = TRUE)) {
  library(TMBhelper)
}

# --- Source Model Functions ---
source(file.path(rtmb_dir, "R", "utilities.R"))
source(file.path(rtmb_dir, "R", "model_funs.R"))
source(file.path(rtmb_dir, "R", "Rpm.R"))

# --- Load ADMB Results for Comparison ---
library(ebswp)

# Path to ADMB rep file
admb_rep_path <- file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.rep")
if (!file.exists(admb_rep_path)) {
  # Try alternative location
  admb_rep_path <- file.path(pollock_root, "runs", "for_rtmb", "pm.rep")
}

if (file.exists(admb_rep_path)) {
  pm <- read_rep(admb_rep_path)
  cat("Loaded ADMB rep file:", admb_rep_path, "\n")

  # Clean up ADMB arrays (indexing differences)
  pm$phat_ats <- pm$phat_ats[, 2:16]
  pm$phat_bts <- pm$phat_bts[, 2:16]
  pm$sel_ats <- pm$sel_ats[31:61, ]
  pm$sel_bts <- pm$sel_bts[19:61, ]
  pm$phat_fsh <- pm$phat_fsh[, 2:16]
  pm$bts_like <- pm$surv_like[1]
  pm$ats_like <- pm$surv_like[2]
  pm$ats_age1_like <- pm$surv_like[3]
  pm$SSB <- pm$SSB[, 2]
} else {
  warning("ADMB rep file not found at: ", admb_rep_path)
  pm <- NULL
}

# --- Load Parameters from ADMB Converged Solution ---
admb_par_path <- file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.par")
if (!file.exists(admb_par_path)) {
  admb_par_path <- file.path(pollock_root, "runs", "rtmb", "pm.par")
}

if (file.exists(admb_par_path)) {
  parms <- read_pars(admb_par_path)
  parms$steepness <- 0.67
  cat("Loaded ADMB par file:", admb_par_path, "\n")
  cat("Set RTMB steepness to fixed value:", parms$steepness, "\n")
} else {
  stop("ADMB par file not found. Required for RTMB initialization.")
}

# --- Load Data ---
data <- Get_Data()
if (!is.null(pm)) {
  data$sam_bts <- floor(pm$sam_bts)  # Avoid non-integer sample sizes
}

# --- Parameter Mapping ---
# Fix parameters that are not being estimated in this configuration

fixed_params <- c(
  # Weight-at-age parameters
  "log_K", "d_scale", "L1", "L2",
  # Initial conditions (use ADMB values)
  "log_avginit", "log_avgrec", "log_avg_F",
  # BTS catchability (estimated analytically)
  "log_q_bts",
  # Future projections
  "rec_dev_future",
  # Natural mortality
  "natmort_phi",
  # Larval transport
  "larv_rec_devs",
  # BTS selectivity (using logistic)
  "sel_devs_bts", "sel_slp_bts", "sel_a50_bts", "sel_age_one_bts",
  # Fishery selectivity (using coefficients, not logistic)
  "sel_trm2_fsh", "sel_dif1_fsh", "sel_a501_fsh", "sel_trm1_fsh",
  "sel_dif2_fsh", "sel_dif1_fsh_dev", "sel_a501_fsh_dev", "sel_trm2_fsh_dev",
  # Temperature effects
  "resid_temp_x1", "resid_temp_x2",
  # Other fixed
  "log_q_std_area", "bt_slope", "sigr", "steepness",
  "sel_coffs_bts",
  # M deviations
  "M_dev",
  # Predation parameters
  "log_a_II", "log_b_II", "log_a_II_vec", "log_b_II_vec",
  "log_rho", "log_resid_M",
  # Weight-at-age
  "log_alpha"
)

map_obj <- create_map_from_par(
  parms, parms,
  exact_names = fixed_params,
  exclude_patterns = "xxx"
)

# --- Create TMB Object ---
obj <- MakeADFun(rpm, parms, map = map_obj)

cat("\nRTMB model object created successfully\n")
cat("Number of parameters:", length(obj$par), "\n")
cat("Initial objective:", obj$fn(), "\n")

# --- Quick Evaluation at ADMB Values ---
if (exists("return_nll_only")) {
  # If running in evaluation mode
  cat("Evaluating at ADMB parameter values...\n")
  rtmb_result <- rpm(parms)
  cat("RTMB NLL:", rtmb_result$tot_like, "\n")
  if (!is.null(pm)) {
    cat("ADMB NLL:", sum(pm$NLL), "\n")
    cat("Difference:", rtmb_result$tot_like - sum(pm$NLL), "\n")
  }
}
