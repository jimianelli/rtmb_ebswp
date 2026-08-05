#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
rtmb_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(rtmb_root, "R", "Rpm.R"))) stop("Run from repository root")
if (!nzchar(Sys.getenv("POLLOCK_ROOT"))) stop("Set POLLOCK_ROOT")
source_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
out_dir <- file.path(rtmb_root, "analysis", "output", "double_logistic_experiments")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_fresh_env <- function() {
  e <- new.env(parent = globalenv())
  e$rm <- function(...) invisible(NULL)
  e$source <- function(file, ...) base::source(file, local = parent.frame(), ...)
  base::source(file.path(rtmb_root, "R", "utils-rtmb.R"), local = e)
  base::source(file.path(rtmb_root, "R", "config.R"), local = e)
  e
}

fixed_names <- function(e) c(
  "log_K", "d_scale", "L1", "L2", "log_avginit", "log_avgrec",
  "log_avg_F", "log_q_bts", "rec_dev_future", "natmort_phi",
  "larv_rec_devs", "sel_devs_bts", "sel_slp_bts", "sel_a50_bts",
  "sel_age_one_bts", "sel_trm2_fsh", "sel_dif1_fsh", "sel_a501_fsh",
  "sel_trm1_fsh", "sel_dif2_fsh", "sel_dif1_fsh_dev",
  "sel_a501_fsh_dev", "sel_trm2_fsh_dev",
  e$inactive_fishery_selectivity_parameters(2L),
  "resid_temp_x1", "resid_temp_x2",
  "log_q_std_area", "bt_slope",
  "sigr", "steepness", "sel_coffs_bts", "sel_tv_ar1_weight_fsh", "M_dev",
  "log_a_II", "log_b_II", "log_a_II_vec", "log_b_II_vec", "log_rho",
  "log_resid_M", "log_alpha"
)

build_model <- function(prior_scale, rw_scale, start_matrix = NULL, pooled = FALSE,
                        random_effects = FALSE, hierarchical_stage = 0L,
                        mean_start = NULL, dev_start = NULL, cv = 0.05,
                        old_age_cap = TRUE) {
  e <- make_fresh_env()
  data <- e$data
  data$fishery_sel_form <- 2L
  data$return_nll_only <- 1L
  data$sel_double_logistic_prior_scale <- prior_scale
  data$sel_double_logistic_rw_scale <- rw_scale
  data$sel_double_logistic_hierarchical <- as.integer(hierarchical_stage > 0L)
  data$sel_double_logistic_cv <- cv
  data$fishery_sel_old_age_cap <- as.integer(old_age_cap)
  if (is.null(data$fishery_sel_spline_basis)) {
    data$fishery_sel_spline_basis <- e$make_fishery_sel_spline_basis(
      as.integer(data$nages), as.integer(data$n_fishery_sel_spline_basis %||% 6L)
    )
  }
  parms <- e$add_fishery_selectivity_parameters(e$parms, data)
  nyrs <- as.integer(data$endyr - data$styr + 1L)
  if (is.null(start_matrix)) {
    start_matrix <- matrix(log(c(1.5, 3, 4)), nyrs, 3, byrow = TRUE)
  }
  stopifnot(identical(dim(start_matrix), c(nyrs, 3L)))
  parms$sel_double_logistic_fsh <- start_matrix
  if (is.null(mean_start)) mean_start <- log(c(1.5, 3, 4))
  parms$sel_double_logistic_mean_fsh <- mean_start
  if (is.null(dev_start)) dev_start <- matrix(0, nyrs, 3)
  stopifnot(identical(dim(dev_start), c(nyrs, 3L)))
  parms$sel_double_logistic_dev_fsh <- dev_start
  fixed <- c(fixed_names(e), "sel_double_logistic_mean_fsh",
             "sel_double_logistic_dev_fsh")
  if (hierarchical_stage > 0L) {
    fixed <- setdiff(fixed, "sel_double_logistic_mean_fsh")
    fixed <- c(fixed, "sel_double_logistic_fsh")
  }
  if (hierarchical_stage == 2L) {
    fixed <- setdiff(fixed, "sel_double_logistic_dev_fsh")
  }
  map_obj <- e$create_map_from_par(parms, parms, exact_names = fixed,
                                   exclude_patterns = "xxx")
  if (pooled) map_obj$sel_double_logistic_fsh <- factor(rep(1:3, each = nyrs))
  e$data <- data
  rpm <- e$rpm
  environment(rpm) <- e
  random <- if (random_effects && hierarchical_stage == 2L) {
    "sel_double_logistic_dev_fsh"
  } else if (random_effects) {
    "sel_double_logistic_fsh"
  } else NULL
  obj <- RTMB::MakeADFun(rpm, parms, map = map_obj, data = data, random = random)
  list(e = e, data = data, rpm = rpm, obj = obj, nyrs = nyrs)
}

make_bounds <- function(obj) {
  lower <- rep(-Inf, length(obj$par)); upper <- rep(Inf, length(obj$par))
  idx <- which(names(obj$par) == "sel_double_logistic_fsh")
  if (length(idx)) {
    n_each <- length(idx) / 3L
    lower[idx] <- rep(log(c(0.25, 1.0, 1.5)), each = n_each)
    upper[idx] <- rep(log(c(5.0, 8.0, 25.0)), each = n_each)
  }
  list(lower = lower, upper = upper)
}

fit_model <- function(name, prior_scale, rw_scale, start_matrix = NULL,
                      pooled = FALSE, seed = NA_integer_, random_effects = FALSE,
                      hierarchical_stage = 0L, mean_start = NULL,
                      dev_start = NULL, cv = 0.05,
                      old_age_cap = TRUE) {
  cat("\n--- ", name, " ---\n", sep = "")
  b <- build_model(
    prior_scale, rw_scale, start_matrix, pooled, random_effects,
    hierarchical_stage, mean_start, dev_start, cv, old_age_cap
  )
  obj <- b$obj; bounds <- make_bounds(obj); t0 <- Sys.time()
  control1 <- list(eval.max = 10000, iter.max = 6000, rel.tol = 1e-10, x.tol = 1e-10)
  control2 <- list(eval.max = 10000, iter.max = 6000, rel.tol = 1e-12, x.tol = 1e-12)
  fit1 <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds$lower,
                 upper = bounds$upper, control = control1)
  fit <- nlminb(fit1$par, obj$fn, obj$gr, lower = bounds$lower,
                upper = bounds$upper, control = control2)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fitted <- obj$env$parList(fit$par)
  b$e$data$return_nll_only <- 0L
  report <- b$rpm(fitted)
  idx <- which(names(obj$par) == "sel_double_logistic_fsh")
  at_bound <- if (length(idx)) sum(
    abs(fit$par[idx] - bounds$lower[idx]) < 1e-5 |
      abs(fit$par[idx] - bounds$upper[idx]) < 1e-5
  ) else 0L
  checkpoint <- list(
    name = name, source_commit = source_commit, seed = seed, pooled = pooled,
    hierarchical_stage = hierarchical_stage, cv = cv,
    old_age_cap = old_age_cap, prior_scale = prior_scale, rw_scale = rw_scale,
    fit_first = fit1, fit = fit, random_effects = random_effects,
    objective = fit$objective, evaluated_total = report$tot_like %||% NA_real_,
    objective_report_gap = if (random_effects) NA_real_ else
      fit$objective - (report$tot_like %||% NA_real_),
    convergence = fit$convergence, message = fit$message,
    max_gradient = max(abs(obj$gr(fit$par)), na.rm = TRUE), at_bound = at_bound,
    hessian_ok = NA, min_hessian_eigen = NA_real_, finite_se = NA,
    seconds = elapsed, fitted_parameters = fitted, report = report
  )
  saveRDS(checkpoint, file.path(out_dir, paste0(name, "_pre_hessian.rds")))
  hess <- tryCatch(optimHess(fit$par, obj$fn, obj$gr), error = identity)
  min_eigen <- NA_real_; hessian_ok <- FALSE; finite_se <- FALSE
  if (is.matrix(hess) && all(is.finite(hess))) {
    eig <- tryCatch(eigen(hess, symmetric = TRUE, only.values = TRUE)$values,
                    error = function(e) NA_real_)
    min_eigen <- suppressWarnings(min(eig, na.rm = TRUE))
    hessian_ok <- is.finite(min_eigen) && min_eigen > 0
    if (hessian_ok) {
      vc <- tryCatch(solve(hess), error = identity)
      finite_se <- is.matrix(vc) && all(is.finite(sqrt(diag(vc))))
    }
  }
  ans <- list(
    name = name, source_commit = source_commit, seed = seed, pooled = pooled,
    hierarchical_stage = hierarchical_stage, cv = cv,
    old_age_cap = old_age_cap,
    prior_scale = prior_scale, rw_scale = rw_scale, fit_first = fit1, fit = fit,
    random_effects = random_effects,
    objective = fit$objective, evaluated_total = report$tot_like %||% NA_real_,
    objective_report_gap = if (random_effects) NA_real_ else
      fit$objective - (report$tot_like %||% NA_real_),
    convergence = fit$convergence, message = fit$message,
    max_gradient = max(abs(obj$gr(fit$par)), na.rm = TRUE), at_bound = at_bound,
    hessian_ok = hessian_ok, min_hessian_eigen = min_eigen,
    finite_se = finite_se, seconds = elapsed,
    fitted_parameters = fitted, report = report
  )
  saveRDS(ans, file.path(out_dir, paste0(name, ".rds")))
  ans
}

as_row <- function(x) tibble(
  scenario = x$name, seed = x$seed, pooled = x$pooled,
  random_effects = x$random_effects %||% FALSE,
  old_age_cap = x$old_age_cap %||% TRUE,
  prior_scale = paste(x$prior_scale, collapse = "/"),
  rw_scale = paste(x$rw_scale, collapse = "/"), objective = x$objective,
  evaluated_total = x$evaluated_total,
  objective_report_gap = x$objective_report_gap,
  convergence = x$convergence, message = x$message,
  max_gradient = x$max_gradient, at_bound = x$at_bound,
  hessian_ok = x$hessian_ok, min_hessian_eigen = x$min_hessian_eigen,
  finite_se = x$finite_se, seconds = x$seconds
)

original_prior <- c(0.75, 0.75, 0.45); original_rw <- c(0.60, 0.60, 0.40)
proposed_prior <- c(0.50, 0.50, 0.35); proposed_rw <- c(0.25, 0.25, 0.20)
mode <- Sys.getenv("FORM2_EXPERIMENT_MODE", "fixed")
results <- list()
if (identical(mode, "random")) {
  metadata <- readRDS(file.path(out_dir, "experiment_metadata.rds"))
  common_start <- metadata$common_start
  for (mult in c(0.50, 0.75, 1.00)) {
    nm <- paste0("random_proposed_", gsub("\\.", "p", format(mult, nsmall = 2)))
    results[[nm]] <- tryCatch(
      fit_model(
        nm, proposed_prior * mult, proposed_rw * mult, common_start,
        random_effects = TRUE
      ),
      error = function(err) list(
        name = nm, seed = NA_integer_, pooled = FALSE, random_effects = TRUE,
        prior_scale = proposed_prior * mult, rw_scale = proposed_rw * mult,
        objective = NA_real_, evaluated_total = NA_real_,
        objective_report_gap = NA_real_, convergence = 99L,
        message = conditionMessage(err), max_gradient = NA_real_,
        at_bound = NA_integer_, hessian_ok = FALSE,
        min_hessian_eigen = NA_real_, finite_se = FALSE, seconds = NA_real_
      )
    )
  }
  summary <- bind_rows(map(results, as_row))
  write_csv(summary, file.path(out_dir, "random_effects_summary.csv"))
} else if (identical(mode, "cv05")) {
  results$cv05_stage1_common <- fit_model(
    "cv05_stage1_common", original_prior, original_rw,
    hierarchical_stage = 1L, cv = 0.05
  )
  mean_hat <- results$cv05_stage1_common$fitted_parameters$sel_double_logistic_mean_fsh
  results$cv05_stage2_random <- tryCatch(
    fit_model(
      "cv05_stage2_random", original_prior, original_rw,
      random_effects = TRUE, hierarchical_stage = 2L,
      mean_start = mean_hat, cv = 0.05
    ),
    error = function(err) list(
      name = "cv05_stage2_random", seed = NA_integer_, pooled = FALSE,
      random_effects = TRUE, hierarchical_stage = 2L, cv = 0.05,
      prior_scale = original_prior, rw_scale = original_rw,
      objective = NA_real_, evaluated_total = NA_real_,
      objective_report_gap = NA_real_, convergence = 99L,
      message = conditionMessage(err), max_gradient = NA_real_,
      at_bound = NA_integer_, hessian_ok = FALSE,
      min_hessian_eigen = NA_real_, finite_se = FALSE, seconds = NA_real_
    )
  )
  summary <- bind_rows(map(results, as_row))
  write_csv(summary, file.path(out_dir, "cv05_summary.csv"))
} else if (identical(mode, "cv05_no_old_age_cap")) {
  results$cv05_no_old_age_cap_stage1_common <- fit_model(
    "cv05_no_old_age_cap_stage1_common", original_prior, original_rw,
    hierarchical_stage = 1L, cv = 0.05, old_age_cap = FALSE
  )
  mean_hat <- results$cv05_no_old_age_cap_stage1_common$fitted_parameters$sel_double_logistic_mean_fsh
  results$cv05_no_old_age_cap_stage2_random <- tryCatch(
    fit_model(
      "cv05_no_old_age_cap_stage2_random", original_prior, original_rw,
      random_effects = TRUE, hierarchical_stage = 2L,
      mean_start = mean_hat, cv = 0.05, old_age_cap = FALSE
    ),
    error = function(err) list(
      name = "cv05_no_old_age_cap_stage2_random", seed = NA_integer_,
      pooled = FALSE, random_effects = TRUE, hierarchical_stage = 2L,
      cv = 0.05, old_age_cap = FALSE,
      prior_scale = original_prior, rw_scale = original_rw,
      objective = NA_real_, evaluated_total = NA_real_,
      objective_report_gap = NA_real_, convergence = 99L,
      message = conditionMessage(err), max_gradient = NA_real_,
      at_bound = NA_integer_, hessian_ok = FALSE,
      min_hessian_eigen = NA_real_, finite_se = FALSE, seconds = NA_real_
    )
  )
  summary <- bind_rows(map(results, as_row))
  write_csv(summary, file.path(out_dir, "cv05_no_old_age_cap_summary.csv"))
} else if (identical(mode, "cv_ladder_no_old_age_cap")) {
  stage1_path <- file.path(out_dir, "cv05_no_old_age_cap_stage1_common.rds")
  if (file.exists(stage1_path)) {
    results$stage1_common <- readRDS(stage1_path)
  } else {
    results$stage1_common <- fit_model(
      "cv_ladder_stage1_common", original_prior, original_rw,
      hierarchical_stage = 1L, cv = 0.05, old_age_cap = FALSE
    )
  }
  mean_hat <- results$stage1_common$fitted_parameters$sel_double_logistic_mean_fsh
  dev_hat <- matrix(0, nrow(results$stage1_common$report$sel_fsh), 3)
  warm_path <- Sys.getenv(
    "FORM2_CV_WARM_START",
    file.path(out_dir, "cv0p075_no_old_age_cap_stage2_random.rds")
  )
  if (file.exists(warm_path)) {
    warm <- readRDS(warm_path)
    mean_hat <- warm$fitted_parameters$sel_double_logistic_mean_fsh
    dev_hat <- warm$fitted_parameters$sel_double_logistic_dev_fsh
  }
  cv_values <- as.numeric(strsplit(
    Sys.getenv("FORM2_CV_LADDER", "0.075,0.10,0.15,0.20"), ",",
    fixed = TRUE
  )[[1]])
  if (!length(cv_values) || any(!is.finite(cv_values)) || any(cv_values <= 0)) {
    stop("FORM2_CV_LADDER must contain positive comma-separated CV values")
  }
  for (cv_value in cv_values) {
    cv_tag <- gsub("\\.", "p", format(cv_value, trim = TRUE, scientific = FALSE))
    nm <- paste0("cv", cv_tag, "_no_old_age_cap_stage2_random")
    results[[nm]] <- tryCatch(
      fit_model(
        nm, original_prior, original_rw,
        random_effects = TRUE, hierarchical_stage = 2L,
        mean_start = mean_hat, dev_start = dev_hat,
        cv = cv_value, old_age_cap = FALSE
      ),
      error = function(err) list(
        name = nm, seed = NA_integer_, pooled = FALSE,
        random_effects = TRUE, hierarchical_stage = 2L,
        cv = cv_value, old_age_cap = FALSE,
        prior_scale = original_prior, rw_scale = original_rw,
        objective = NA_real_, evaluated_total = NA_real_,
        objective_report_gap = NA_real_, convergence = 99L,
        message = conditionMessage(err), max_gradient = NA_real_,
        at_bound = NA_integer_, hessian_ok = FALSE,
        min_hessian_eigen = NA_real_, finite_se = FALSE, seconds = NA_real_
      )
    )
    if (!is.null(results[[nm]]$fitted_parameters$sel_double_logistic_mean_fsh)) {
      mean_hat <- results[[nm]]$fitted_parameters$sel_double_logistic_mean_fsh
      dev_hat <- results[[nm]]$fitted_parameters$sel_double_logistic_dev_fsh
    }
  }
  summary <- bind_rows(map(results[names(results) != "stage1_common"], as_row))
  write_csv(summary, file.path(out_dir, "cv_ladder_no_old_age_cap_summary.csv"))
} else {
  results$original_reproduction <- fit_model("original_reproduction", original_prior, original_rw)
  results$pooled_original <- fit_model("pooled_original", original_prior, original_rw, pooled = TRUE)
  pooled_triplet <- results$pooled_original$fitted_parameters$sel_double_logistic_fsh[1, ]
  nyrs <- nrow(results$pooled_original$fitted_parameters$sel_double_logistic_fsh)
  common_start <- matrix(pooled_triplet, nyrs, 3, byrow = TRUE)
  results$common_start_original <- fit_model("common_start_original", original_prior, original_rw, common_start)
  results$proposed_prior <- fit_model("proposed_prior", proposed_prior, proposed_rw, common_start)
  for (mult in c(0.75, 1.25, 1.50)) {
    nm <- paste0("prior_multiplier_", gsub("\\.", "p", format(mult, nsmall = 2)))
    results[[nm]] <- fit_model(nm, proposed_prior * mult, proposed_rw * mult, common_start)
  }
  for (prior_name in c("original", "proposed")) {
    ps <- if (prior_name == "original") original_prior else proposed_prior
    rs <- if (prior_name == "original") original_rw else proposed_rw
    for (seed in 1:10) {
      set.seed(seed)
      jittered <- common_start + matrix(rnorm(length(common_start), sd = 0.05), nyrs, 3)
      nm <- sprintf("multistart_%s_%02d", prior_name, seed)
      results[[nm]] <- fit_model(nm, ps, rs, jittered, seed = seed)
    }
  }
  summary <- bind_rows(map(results, as_row))
  write_csv(summary, file.path(out_dir, "summary.csv"))
  saveRDS(list(source_commit = source_commit, common_start = common_start,
               summary = summary), file.path(out_dir, "experiment_metadata.rds"))
}
print(summary, n = Inf)
cat("\nResults written to ", out_dir, "\n", sep = "")
