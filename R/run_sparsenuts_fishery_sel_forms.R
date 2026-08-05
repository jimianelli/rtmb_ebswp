#!/usr/bin/env Rscript

# Run SparseNUTS MCMC for fishery selectivity scenarios. Form 2 uses the
# accepted two-stage hierarchical double-logistic configuration.
# and write plots (marginals + pairs) for slow-mixing parameters.
#
# Usage:
#   Rscript R/run_sparsenuts_fishery_sel_forms.R
#   POLLOCK_ROOT=/path/to/pollock Rscript R/run_sparsenuts_fishery_sel_forms.R

suppressPackageStartupMessages({
  library(RTMB)
  library(tidyverse)
  library(SparseNUTS)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

rtmb_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(rtmb_dir, "R", "config.R"))) {
  stop("Run from repo root (directory containing R/config.R)")
}

resolve_pollock_root <- function(rtmb_dir) {
  has_bridge <- function(path) {
    if (is.na(path) || !nzchar(path)) return(FALSE)
    file.exists(file.path(path, "admb", "runs", "for_rtmb", "pm.rep")) &&
      file.exists(file.path(path, "admb", "runs", "for_rtmb", "pm.par")) &&
      file.exists(file.path(path, "admb", "runs", "data", "pm_24.dat"))
  }
  env_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
  if (is.na(env_root) || !nzchar(env_root)) {
    env_root <- Sys.getenv("POLLOCK_BASE", unset = NA_character_)
  }
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }
  candidates <- c(rtmb_dir, dirname(rtmb_dir), file.path(dirname(rtmb_dir), "pollock"))
  for (cand in candidates) {
    if (has_bridge(cand)) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  stop("Cannot locate pollock bridge inputs. Set POLLOCK_ROOT or use bundled admb/runs files.")
}
pollock_root <- resolve_pollock_root(rtmb_dir)

forms <- as.integer(strsplit(Sys.getenv("SPARSENUTS_FORMS", "2"), ",", fixed = TRUE)[[1]])
chains <- as.integer(Sys.getenv("SPARSENUTS_CHAINS", "8"))
cores <- as.integer(Sys.getenv("SPARSENUTS_CORES", as.character(min(chains, 8L))))
iter <- as.integer(Sys.getenv("SPARSENUTS_ITER", "2000"))
warmup <- as.integer(Sys.getenv("SPARSENUTS_WARMUP", as.character(floor(iter / 2))))
seed <- as.integer(Sys.getenv("SPARSENUTS_SEED", "123"))
form2_cv <- as.numeric(Sys.getenv("FORM2_CV", "0.20"))
build_only <- tolower(Sys.getenv("SPARSENUTS_BUILD_ONLY", "false")) %in%
  c("1", "true", "yes")
stopifnot(length(forms) > 0L, all(forms %in% c(0L, 2L, 5L)))
stopifnot(chains > 0L, cores > 0L, iter > 1L, warmup >= 0L, warmup < iter)
out_dir <- file.path("analysis", "output", "sparsenuts", "fishery_sel_forms")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Build an RTMB obj for a given fishery_sel_form by sourcing config.R into an env
make_obj_for_form <- function(form) {
  rtmb_env <- new.env(parent = globalenv())
  rtmb_env$rm <- function(...) invisible(NULL)
  rtmb_env$source <- function(file, ...) base::source(file, local = parent.frame(), ...)

  # Source config to get data/parms/rpm and helper fns
  source(file.path(rtmb_dir, "R", "config.R"), local = rtmb_env)
  # Pull only the fishery-selectivity helper functions from utils-rtmb.R
  # (sourcing the whole file can overwrite other model functions used in rpm()).
  tmp_env <- new.env(parent = baseenv())
  base::source(file.path(rtmb_dir, "R", "utils-rtmb.R"), local = tmp_env)
  rtmb_env$add_fishery_selectivity_parameters <- tmp_env$add_fishery_selectivity_parameters
  rtmb_env$inactive_fishery_selectivity_parameters <- tmp_env$inactive_fishery_selectivity_parameters

  # Override form in data (used by inactive_fishery_selectivity_parameters())
  rtmb_env$data$fishery_sel_form <- as.integer(form)

  if (form == 2L) {
    rtmb_env$data$sel_double_logistic_hierarchical <- 1L
    rtmb_env$data$sel_double_logistic_cv <- form2_cv
    rtmb_env$data$fishery_sel_old_age_cap <- 0L
  }

  # Ensure parameter blocks exist
  rtmb_env$parms <- rtmb_env$add_fishery_selectivity_parameters(rtmb_env$parms, rtmb_env$data)

  start_file <- Sys.getenv(
    "FORM2_START_FILE",
    file.path(rtmb_dir, "analysis", "output", "double_logistic_experiments",
              "cv0p2_no_old_age_cap_stage2_random.rds")
  )
  if (form == 2L && file.exists(start_file)) {
    saved <- readRDS(start_file)
    fitted <- saved$fitted_parameters
    for (nm in intersect(names(fitted), names(rtmb_env$parms))) {
      if (identical(dim(fitted[[nm]]), dim(rtmb_env$parms[[nm]])) &&
          length(fitted[[nm]]) == length(rtmb_env$parms[[nm]])) {
        rtmb_env$parms[[nm]] <- fitted[[nm]]
      }
    }
  }

  # Recreate map_obj + obj exactly as config.R does, but respecting new form
  fixed_params <- c(
    "log_K", "d_scale", "L1", "L2",
    "log_avginit", "log_avgrec", "log_avg_F",
    "log_q_bts",
    "rec_dev_future",
    "natmort_phi",
    "larv_rec_devs",
    "sel_devs_bts", "sel_slp_bts", "sel_a50_bts", "sel_age_one_bts",
    "sel_trm2_fsh", "sel_dif1_fsh", "sel_a501_fsh", "sel_trm1_fsh",
    "sel_dif2_fsh", "sel_dif1_fsh_dev", "sel_a501_fsh_dev", "sel_trm2_fsh_dev",
    rtmb_env$inactive_fishery_selectivity_parameters(rtmb_env$data$fishery_sel_form),
    "resid_temp_x1", "resid_temp_x2",
    "log_q_std_area", "bt_slope", "sigr", "steepness",
    "sel_coffs_bts",
    if (form == 2L) "sel_double_logistic_fsh" else character(0),
    "M_dev",
    "log_a_II", "log_b_II", "log_a_II_vec", "log_b_II_vec",
    "log_rho", "log_resid_M",
    "log_alpha"
  )

  map_obj <- rtmb_env$create_map_from_par(
    rtmb_env$parms, rtmb_env$parms,
    exact_names = fixed_params,
    exclude_patterns = "xxx"
  )

  random <- if (form == 2L) "sel_double_logistic_dev_fsh" else NULL
  obj <- MakeADFun(rtmb_env$rpm, rtmb_env$parms, map = map_obj, random = random)

  list(obj = obj, data = rtmb_env$data, parms = rtmb_env$parms)
}

run_one <- function(form) {
  cat("\n=== SparseNUTS: fishery_sel_form=", form, " ===\n", sep = "")
  built <- make_obj_for_form(form)
  obj <- built$obj
  data <- built$data

  if (build_only) {
    objective <- obj$fn(obj$par)
    gradient <- obj$gr(obj$par)
    cat("Build-only check passed; objective =", objective,
        "; max fixed-effect gradient =", max(abs(gradient)), "\n")
    return(invisible(list(objective = objective, max_gradient = max(abs(gradient)))))
  }

  fit_file <- file.path(out_dir, sprintf("rtmb_ebswp_sparsenuts_form_%d.rds", form))

  # Run MCMC
  cat("Configuration:", chains, "chains;", iter, "iterations;", warmup,
      "warmup;", cores, "cores.\n")
  fit <- SparseNUTS::sample_snuts(
    obj,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    seed = seed,
    init = "last.par.best",
    metric = "diag",
    globals = list(data = data)
  )
  attr(fit, "rtmb_ebswp_sparsenuts") <- list(
    fishery_sel_form = form,
    hierarchical_form2 = form == 2L,
    form2_cv = if (form == 2L) form2_cv else NA_real_,
    old_age_cap = if (form == 2L) FALSE else NA,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    created = Sys.time(),
    execution = "serial chains requested with cores=1"
  )
  saveRDS(fit, fit_file)
  cat("Saved fit: ", fit_file, "\n", sep = "")

  # Identify slow parameters from monitor
  slow_idx <- integer(0)
  if (!is.null(fit$monitor) && "ess_bulk" %in% names(fit$monitor)) {
    ess <- fit$monitor$ess_bulk
    slow_idx <- order(ess)[seq_len(min(6, length(ess)))]
  }

  # Plot pairs for slow params (SparseNUTS S3 method)
  pairs_png <- file.path(fig_dir, sprintf("sparsenuts_pairs_slow_form_%d.png", form))
  png(pairs_png, width = 1800, height = 1600, res = 180)
  try({
    if (length(slow_idx) >= 2) {
      pairs(fit, pars = slow_idx, order = "orig", diag = "trace")
    } else {
      plot.new(); text(0.5, 0.5, "No slow-parameter indices available")
    }
  }, silent = TRUE)
  dev.off()

  # Also generate pairs.tmbfit explicitly (same underlying function, for your request)
  pairs_tmbfit_png <- file.path(fig_dir, sprintf("sparsenuts_pairs_tmbfit_form_%d.png", form))
  png(pairs_tmbfit_png, width = 1800, height = 1600, res = 180)
  try({
    if (length(slow_idx) >= 2) {
      SparseNUTS:::pairs.tmbfit(fit, pars = slow_idx, order = "orig", diag = "trace")
    } else {
      plot.new(); text(0.5, 0.5, "No slow-parameter indices available")
    }
  }, silent = TRUE)
  dev.off()

  # Marginals for slow params
  marg_png <- file.path(fig_dir, sprintf("sparsenuts_marginals_slow_form_%d.png", form))
  png(marg_png, width = 1800, height = 1400, res = 180)
  try({
    if (length(slow_idx) >= 1) {
      SparseNUTS::plot_marginals(fit, pars = slow_idx, nrow = 3)
    } else {
      plot.new(); text(0.5, 0.5, "No slow-parameter indices available")
    }
  }, silent = TRUE)
  dev.off()

  invisible(list(fit_file = fit_file, pairs_png = pairs_png, marg_png = marg_png, pairs_tmbfit_png = pairs_tmbfit_png))
}

results <- lapply(forms, run_one)
cat("\nDone. Outputs in: ", out_dir, "\n", sep = "")
