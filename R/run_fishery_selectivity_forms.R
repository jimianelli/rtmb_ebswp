#!/usr/bin/env Rscript

# Run alternative fishery selectivity forms (real model fits) and compare to base.
#
# Usage:
#   Rscript R/run_fishery_selectivity_forms.R
#   POLLOCK_ROOT=~/workspace/pollock Rscript R/run_fishery_selectivity_forms.R
#
# Outputs:
#   analysis/output/fishery_sel_forms/fishery_sel_form_<id>.rds
#   analysis/output/fishery_sel_forms/summary.csv

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

rtmb_root <- normalizePath(getwd(), mustWork = TRUE)

output_root <- Sys.getenv(
  "SELECTIVITY_OUTPUT_ROOT",
  file.path("analysis", "output")
)
out_dir <- file.path(rtmb_root, output_root, "fishery_sel_forms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
form5_random <- tolower(Sys.getenv("FORM5_RANDOM_EFFECTS", "false")) %in%
  c("1", "true", "yes")
form5_estimate_hyper <- tolower(Sys.getenv("FORM5_ESTIMATE_HYPER", "false")) %in%
  c("1", "true", "yes")
form5_fix_sigma <- tolower(Sys.getenv("FORM5_FIX_SIGMA", "false")) %in%
  c("1", "true", "yes")
form5_random_field_only <- tolower(Sys.getenv(
  "FORM5_RANDOM_FIELD_ONLY", "false"
)) %in% c("1", "true", "yes")
form5_max_age <- as.integer(Sys.getenv("FORM5_MAX_AGE", "15"))

# Build a completely fresh model environment for each form to avoid any RTMB/env contamination.
make_fresh_env <- function() {
  e <- new.env(parent = globalenv())
  e$rm <- function(...) invisible(NULL)
  e$source <- function(file, ...) base::source(file, local = parent.frame(), ...)
  base::source(file.path(rtmb_root, "R", "utils-rtmb.R"), local = e)
  base::source(file.path(rtmb_root, "R", "config.R"), local = e)
  e
}

fit_one <- function(form_id, label) {
  form_id <- as.integer(form_id)
  cat("\n--- Fitting fishery_sel_form=", form_id, " (", label, ") ---\n", sep = "")

  e <- make_fresh_env()

  # Overwrite the selectivity form flag and rebuild the objective so mapping reflects it.
  data <- e$data
  data$fishery_sel_form <- form_id
  data$return_nll_only <- 1
  if (is.null(data$fishery_sel_spline_basis)) {
    data$fishery_sel_spline_basis <- e$make_fishery_sel_spline_basis(
      nages = as.integer(data$nages),
      nbasis = as.integer(data$n_fishery_sel_spline_basis %||% 6L)
    )
  }

  parms <- e$add_fishery_selectivity_parameters(e$parms, data)

  # Nudge starting values for alternative selectivity forms so the optimizer
  # doesn't get stuck at the zero / no-time-variation mode.
  if (form_id == 2L) {
    # 3-parameter double-logistic: one p1/p2/p3 triplet per model year.
    nyrs <- as.integer(data$endyr - data$styr + 1)
    parms$sel_double_logistic_fsh <- matrix(
      log(c(1.5, 3, 4)),
      nrow = nyrs,
      ncol = 3,
      byrow = TRUE
    )
  }
  if (form_id == 5L) {
    set.seed(123)
    # Larger starting sigma encourages non-zero time-varying field
    parms$log_sel_tv_ar1_sigma_fsh <- log(1.2)
    # Relax the AR1 penalty weight (default is 1.0 inside the model).
    # This is model data, not an estimated parameter.
    data$sel_tv_ar1_weight_fsh <- 0.25
    # Keep rho near 0 initially (on working scale)
    parms$sel_tv_ar1_rho_fsh <- c(0, 0)
    # Seed a small non-zero AR1 field (year x age)
    nyrs <- as.integer(data$endyr - data$styr + 1)
    nages <- min(as.integer(data$nages), form5_max_age)
    if (nages < 2L) stop("FORM5_MAX_AGE must be at least 2")
    data$fishery_sel_old_age_cap <- as.integer(nages < as.integer(data$nages))
    parms$sel_tv_ar1_fsh <- matrix(rnorm(nyrs * nages, sd = 0.05), nrow = nyrs, ncol = nages)
    if (form5_random) {
      data$fishery_sel_old_age_cap <- 0L
      data$sel_tv_ar1_weight_fsh <- 1
      fixed_start <- file.path(out_dir, "fishery_sel_form_5.rds")
      if (file.exists(fixed_start)) {
        saved <- readRDS(fixed_start)$fixed_parameters
        compatible <- intersect(names(saved), names(parms))
        for (nm in compatible) {
          same_length <- length(saved[[nm]]) == length(parms[[nm]])
          same_dims <- identical(dim(saved[[nm]]), dim(parms[[nm]]))
          if (same_length && same_dims) parms[[nm]] <- saved[[nm]]
        }
      }
      if (form5_fix_sigma) parms$log_sel_tv_ar1_sigma_fsh <- log(0.2368)
    }
  }

  # rpm() reads `data` from its function environment via getAll(); make the
  # scenario data visible before taping, optimization, and reporting.
  e$data <- data

  fixed_params <- c(
    # Weight-at-age parameters
    "log_K", "d_scale", "L1", "L2",
    # Initial conditions
    "log_avginit", "log_avgrec", "log_avg_F",
    # BTS catchability
    "log_q_bts",
    # Future projections
    "rec_dev_future",
    # Natural mortality
    "natmort_phi",
    # Larval transport
    "larv_rec_devs",
    # BTS selectivity (inactive base terms)
    "sel_devs_bts", "sel_slp_bts", "sel_a50_bts", "sel_age_one_bts",
    # Fishery logistic terms (still dormant)
    "sel_trm2_fsh", "sel_dif1_fsh", "sel_a501_fsh", "sel_trm1_fsh",
    "sel_dif2_fsh", "sel_dif1_fsh_dev", "sel_a501_fsh_dev", "sel_trm2_fsh_dev",
    # Fix all fishery selectivity parameter blocks that are inactive under this form
    e$inactive_fishery_selectivity_parameters(data$fishery_sel_form),
    if (form_id == 5L && form5_random && !form5_estimate_hyper) {
      c("sel_tv_ar1_rho_fsh", "log_sel_tv_ar1_sigma_fsh")
    } else if (form_id == 5L && form5_random && form5_fix_sigma) {
      "log_sel_tv_ar1_sigma_fsh"
    } else character(0),
    # Other fixed
    "log_q_std_area", "bt_slope", "sigr", "steepness",
    "sel_coffs_bts",
    "sel_tv_ar1_weight_fsh",
    "M_dev",
    "log_a_II", "log_b_II", "log_a_II_vec", "log_b_II_vec",
    "log_rho", "log_resid_M",
    "log_alpha"
  )
  if (form_id == 5L && form5_random && form5_random_field_only) {
    # Diagnostic checkpoint: integrate only the AR1 field while holding every
    # outer/model parameter at the saved fixed-effects Form 5 estimates.
    fixed_params <- setdiff(names(parms), "sel_tv_ar1_fsh")
  }

  map_obj <- e$create_map_from_par(
    parms, parms,
    exact_names = fixed_params,
    exclude_patterns = "xxx"
  )

  rpm <- e$rpm
  environment(rpm) <- e
  random <- if (form_id == 5L && form5_random) "sel_tv_ar1_fsh" else NULL
  obj <- RTMB::MakeADFun(rpm, parms, map = map_obj, data = data, random = random)

  lower <- rep(-Inf, length(obj$par))
  upper <- rep(Inf, length(obj$par))
  if (form_id == 2L) {
    idx <- which(names(obj$par) == "sel_double_logistic_fsh")
    if (length(idx) > 0) {
      nyrs <- length(idx) / 3
      lower[idx] <- rep(log(c(0.25, 1.0, 1.5)), each = nyrs)
      upper[idx] <- rep(log(c(5.0, 8.0, 25.0)), each = nyrs)
    }
  }

  t0 <- Sys.time()
  if (length(obj$par) == 0L) {
    fit <- list(
      par = obj$par,
      objective = obj$fn(obj$par),
      convergence = 0L,
      message = "Random-field diagnostic with all outer parameters fixed"
    )
  } else {
    fit <- nlminb(
      obj$par, obj$fn, obj$gr,
      lower = lower,
      upper = upper,
      control = list(eval.max = 5000, iter.max = 3000)
    )
  }
  t1 <- Sys.time()

  # parList() expects the outer/fixed vector. For random-effects models it
  # inserts that vector into the full parameter vector and retains the latest
  # conditional modes for the random block.
  if (!is.null(random)) obj$fn(fit$par)
  fitted_parms <- obj$env$parList(fit$par)
  e$data$return_nll_only <- 0
  report <- rpm(fitted_parms)
  common_input_names <- c(
    "oac_fsh", "oac_bts", "oac_ats", "sam_fsh", "sam_bts", "sam_ats",
    "yrs_fsh_data", "yrs_bts_data", "yrs_ats_data"
  )
  for (nm in common_input_names) report[[nm]] <- data[[nm]]
  evaluated_total <- report$tot_like %||% NA_real_

  res <- list(
    form = form_id,
    label = label,
    bridge_case = Sys.getenv("EBSWP_BRIDGE_CASE", "legacy_pm_bridge"),
    admb_run_dir = Sys.getenv("EBSWP_ADMB_RUN_DIR", "runs/rtmb"),
    bts_comp_normalization = Sys.getenv(
      "EBSWP_BTS_COMP_NORMALIZATION", "legacy_age2plus"
    ),
    objective = fit$objective,
    evaluated_total = evaluated_total,
    convergence = fit$convergence,
    message = fit$message,
    n_parameters = length(fit$par),
    max_gradient = if (length(fit$par)) {
      max(abs(obj$gr(fit$par)), na.rm = TRUE)
    } else 0,
    seconds = as.numeric(difftime(t1, t0, units = "secs")),
    random_effects = !is.null(random),
    random_effect_names = random,
    old_age_cap = if (form_id == 5L) data$fishery_sel_old_age_cap %||% 1L else NA,
    ar1_max_age = if (form_id == 5L) ncol(fitted_parms$sel_tv_ar1_fsh) else NA_integer_,
    ar1_weight = if (form_id == 5L) data$sel_tv_ar1_weight_fsh else NA_real_,
    estimate_ar1_hyperparameters = if (form_id == 5L) {
      form5_estimate_hyper
    } else NA,
    sigma_fixed = if (form_id == 5L) form5_fix_sigma else NA,
    random_field_only = form_id == 5L && form5_random_field_only,
    fixed_parameters = fitted_parms,
    optimizer_parameters = fit$par,
    parameter_table = data.frame(
      term = make.unique(names(fit$par)),
      estimate = as.numeric(fit$par),
      std.error = rep(NA_real_, length(fit$par)),
      gradient = if (length(fit$par)) as.numeric(obj$gr(fit$par)) else numeric(),
      estimation_type = rep("fixed_effect", length(fit$par)),
      module = rep(NA_character_, length(fit$par)),
      year = rep(NA_integer_, length(fit$par)),
      age = rep(NA_integer_, length(fit$par)),
      engine = rep("custom RTMB", length(fit$par)),
      model = rep(label, length(fit$par)),
      stringsAsFactors = FALSE
    ),
    report = report
  )

  suffix <- if (form_id == 5L && form5_max_age < 15L) {
    paste0("_age", form5_max_age, "_plus")
  } else if (form_id == 5L && form5_random) {
    if (form5_random_field_only) {
      "_random_field_only"
    } else if (form5_estimate_hyper && form5_fix_sigma) {
      "_random_sigma_fixed"
    } else if (form5_estimate_hyper) "_random" else "_random_fixed_hyper"
  } else ""
  saveRDS(res, file.path(
    out_dir, sprintf("fishery_sel_form_%d%s.rds", form_id, suffix)
  ))
  res
}

runs <- list(
  list(id = 0L, label = "Base (coefficients + change-year deviations)"),
  list(id = 2L, label = "Time-varying double logistic (annual p1/p2/p3)"),
  list(id = 5L, label = "2D AR1 year×age (all years)" )
)

requested_forms <- as.integer(strsplit(
  Sys.getenv("FISHERY_SEL_FORMS", if (form5_random) "5" else "0,2,5"),
  ",", fixed = TRUE
)[[1]])
runs <- keep(runs, ~ .x$id %in% requested_forms)

results <- purrr::map(runs, ~fit_one(.x$id, .x$label))

summary_tbl <- tibble(
  fishery_sel_form = purrr::map_int(results, "form"),
  label = purrr::map_chr(results, "label"),
  objective = purrr::map_dbl(results, "objective"),
  evaluated_total = purrr::map_dbl(results, "evaluated_total"),
  n_parameters = purrr::map_int(results, "n_parameters"),
  convergence = purrr::map_int(results, "convergence"),
  max_gradient = purrr::map_dbl(results, "max_gradient"),
  seconds = purrr::map_dbl(results, "seconds")
) |>
  arrange(fishery_sel_form)

print(summary_tbl)
readr::write_csv(summary_tbl, file.path(out_dir, "summary.csv"))

cat("\nWrote results to:", out_dir, "\n")
