#!/usr/bin/env Rscript

# Run alternative fishery selectivity forms (real model fits) and compare to base.
#
# Usage:
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

if (is.na(Sys.getenv("POLLOCK_ROOT", unset = NA_character_)) || Sys.getenv("POLLOCK_ROOT") == "") {
  stop("Set POLLOCK_ROOT before running (e.g., POLLOCK_ROOT=~/workspace/pollock)")
}

out_dir <- file.path(rtmb_root, "analysis", "output", "fishery_sel_forms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
    nages <- as.integer(data$nages)
    parms$sel_tv_ar1_fsh <- matrix(rnorm(nyrs * nages, sd = 0.05), nrow = nyrs, ncol = nages)
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
    # Other fixed
    "log_q_std_area", "bt_slope", "sigr", "steepness",
    "sel_coffs_bts",
    "sel_tv_ar1_weight_fsh",
    "M_dev",
    "log_a_II", "log_b_II", "log_a_II_vec", "log_b_II_vec",
    "log_rho", "log_resid_M",
    "log_alpha"
  )

  map_obj <- e$create_map_from_par(
    parms, parms,
    exact_names = fixed_params,
    exclude_patterns = "xxx"
  )

  rpm <- e$rpm
  environment(rpm) <- e
  obj <- RTMB::MakeADFun(rpm, parms, map = map_obj, data = data)

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
  fit <- nlminb(
    obj$par, obj$fn, obj$gr,
    lower = lower,
    upper = upper,
    control = list(eval.max = 5000, iter.max = 3000)
  )
  t1 <- Sys.time()

  fitted_parms <- obj$env$parList(fit$par)
  e$data$return_nll_only <- 0
  report <- rpm(fitted_parms)
  evaluated_total <- report$tot_like %||% NA_real_

  res <- list(
    form = form_id,
    label = label,
    objective = fit$objective,
    evaluated_total = evaluated_total,
    convergence = fit$convergence,
    message = fit$message,
    max_gradient = max(abs(obj$gr(fit$par)), na.rm = TRUE),
    seconds = as.numeric(difftime(t1, t0, units = "secs")),
    fixed_parameters = fitted_parms,
    report = report
  )

  saveRDS(res, file.path(out_dir, sprintf("fishery_sel_form_%d.rds", form_id)))
  res
}

runs <- list(
  list(id = 0L, label = "Base (coefficients + change-year deviations)"),
  list(id = 2L, label = "Time-varying double logistic (annual p1/p2/p3)"),
  list(id = 5L, label = "2D AR1 year×age (all years)" )
)

results <- purrr::map(runs, ~fit_one(.x$id, .x$label))

summary_tbl <- tibble(
  fishery_sel_form = purrr::map_int(results, "form"),
  label = purrr::map_chr(results, "label"),
  objective = purrr::map_dbl(results, "objective"),
  evaluated_total = purrr::map_dbl(results, "evaluated_total"),
  convergence = purrr::map_int(results, "convergence"),
  max_gradient = purrr::map_dbl(results, "max_gradient"),
  seconds = purrr::map_dbl(results, "seconds")
) |>
  arrange(fishery_sel_form)

print(summary_tbl)
readr::write_csv(summary_tbl, file.path(out_dir, "summary.csv"))

cat("\nWrote results to:", out_dir, "\n")
