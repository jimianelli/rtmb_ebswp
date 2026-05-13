#!/usr/bin/env Rscript

# Simulate pseudo-data sets from the saved base-model point estimate.
#
# Usage:
#   NSIM=100 SEED=123 Rscript R/simulate_base_datasets.R
#
# Output:
#   analysis/output/simulated/base_simulated_datasets.rds

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

positive <- function(x, eps = 1e-8) {
  pmax(as.numeric(x), eps)
}

normalize_rows <- function(x, cols = seq_len(ncol(x)), eps = 1e-12) {
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  z <- as.matrix(x[, cols, drop = FALSE])
  z[!is.finite(z) | z < 0] <- 0
  rs <- rowSums(z)
  bad <- !is.finite(rs) | rs <= eps
  if (any(bad)) {
    z[bad, ] <- 1 / length(cols)
    rs[bad] <- 1
  }
  out[, cols] <- sweep(z, 1, rs, "/")
  out
}

rmultinom_rows <- function(prob, size) {
  prob <- as.matrix(prob)
  size <- pmax(0L, as.integer(round(size)))
  out <- matrix(0, nrow = nrow(prob), ncol = ncol(prob), dimnames = dimnames(prob))
  for (i in seq_len(nrow(prob))) {
    out[i, ] <- as.vector(stats::rmultinom(1, size[i], prob[i, ]))
  }
  out
}

rmvnorm_eigen <- function(mu, sigma) {
  mu <- as.numeric(mu)
  sigma <- as.matrix(sigma)
  if (!all(dim(sigma) == length(mu))) {
    stop("Covariance matrix dimensions do not match the mean vector.")
  }
  decomp <- eigen(sigma, symmetric = TRUE)
  values <- pmax(decomp$values, 0)
  z <- stats::rnorm(length(mu))
  mu + as.vector(decomp$vectors %*% (sqrt(values) * z))
}

simulate_lognormal_obs <- function(mu, log_var, eps = 0.01) {
  mu <- positive(mu + eps)
  ans <- exp(stats::rnorm(length(mu), log(mu), sqrt(as.numeric(log_var)))) - eps
  positive(ans)
}

simulate_rtmb_dataset <- function(report, data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  sim <- data
  nages <- as.integer(data$nages)
  mina_ats <- as.integer(data$mina_ats %||% 2L)
  ats_cols <- mina_ats:nages

  catch_sigma <- sqrt(1 / (2 * as.numeric(data$catBio %||% 200)))
  sim$obs_catch <- positive(exp(stats::rnorm(
    length(report$pred_catch),
    log(positive(report$pred_catch)),
    catch_sigma
  )))

  if (!is.null(data$cov_matrix)) {
    sim$ob_bts <- positive(rmvnorm_eigen(report$eb_bts, data$cov_matrix))
  } else {
    sim$ob_bts <- positive(stats::rnorm(
      length(report$eb_bts),
      report$eb_bts,
      as.numeric(data$ob_bts_std)
    ))
  }

  sim$ob_ats <- simulate_lognormal_obs(report$eb_ats, data$lvarb_ats)

  sim$obs_cpue <- positive(stats::rnorm(
    length(report$pred_cpue),
    report$pred_cpue,
    sqrt(as.numeric(data$obs_cpue_var))
  ))

  sim$ob_avo <- positive(stats::rnorm(
    length(report$pred_avo),
    report$pred_avo,
    sqrt(as.numeric(data$obs_avo_var))
  ))

  p_fsh <- normalize_rows(report$phat_fsh)
  p_bts <- normalize_rows(report$phat_bts)
  p_ats <- normalize_rows(report$phat_ats, cols = ats_cols)

  sim$oac_fsh_data <- rmultinom_rows(p_fsh, report$sam_fsh)
  sim$oac_bts <- rmultinom_rows(p_bts, report$sam_bts)
  sim$oac_ats <- rmultinom_rows(p_ats, report$sam_ats)

  sim$oac_fsh <- normalize_rows(sim$oac_fsh_data)
  sim$oac_bts <- normalize_rows(sim$oac_bts)
  sim$oac_ats <- normalize_rows(sim$oac_ats, cols = ats_cols)
  sim$oa1_ats <- simulate_lognormal_obs(
    report$phat_ats[, 1],
    rep((data$age1_sigma_ats %||% 1)^2, nrow(report$phat_ats))
  )
  sim$oac_ats[, 1] <- sim$oa1_ats
  sim$ot_bts <- rowSums(sim$oac_bts[, as.integer(data$mina_bts %||% 2L):nages, drop = FALSE])
  sim$ot_ats <- rowSums(sim$oac_ats[, ats_cols, drop = FALSE])

  attr(sim, "simulation") <- list(
    model = "rtmb_ebswp",
    source = "base point estimate",
    created = Sys.time()
  )
  sim
}

prepare_simulated_estimation_data <- function(sim) {
  sim$return_nll_only <- 1
  sim
}

validate_simulated_estimation_data <- function(sim, report = NULL) {
  required <- c(
    "obs_catch", "ob_bts", "ob_ats", "obs_cpue", "ob_avo",
    "oac_fsh", "oac_bts", "oac_ats",
    "sam_fsh", "sam_bts", "sam_ats",
    "yrs_fsh_data", "yrs_bts_data", "yrs_ats_data",
    "nages", "mina_bts", "mina_ats"
  )
  missing <- setdiff(required, names(sim))
  if (length(missing) > 0) {
    stop("Simulated data are missing required fields: ", paste(missing, collapse = ", "))
  }

  nages <- as.integer(sim$nages)
  if (ncol(sim$oac_fsh) != nages || ncol(sim$oac_bts) != nages || ncol(sim$oac_ats) != nages) {
    stop("Age-composition matrices do not all have nages columns.")
  }

  checks <- c(
    length(sim$obs_catch) == length(sim$styr:sim$endyr),
    length(sim$ob_bts) == length(sim$yrs_bts_data),
    length(sim$ob_ats) == length(sim$yrs_ats_data),
    length(sim$obs_cpue) == length(sim$yrs_cpue),
    length(sim$ob_avo) == length(sim$yrs_avo),
    nrow(sim$oac_fsh) == length(sim$yrs_fsh_data),
    nrow(sim$oac_bts) == length(sim$yrs_bts_data),
    nrow(sim$oac_ats) == length(sim$yrs_ats_data),
    length(sim$sam_fsh) == nrow(sim$oac_fsh),
    length(sim$sam_bts) == nrow(sim$oac_bts),
    length(sim$sam_ats) == nrow(sim$oac_ats)
  )
  if (!all(checks)) {
    stop("Simulated data dimensions are not aligned with year/sample-size fields.")
  }

  row_tol <- 1e-8
  ats_cols <- as.integer(sim$mina_ats):nages
  rows_ok <- c(
    all(abs(rowSums(sim$oac_fsh) - 1) < row_tol),
    all(abs(rowSums(sim$oac_bts) - 1) < row_tol),
    all(abs(rowSums(sim$oac_ats[, ats_cols, drop = FALSE]) - 1) < row_tol)
  )
  if (!all(rows_ok)) {
    stop("Simulated age-composition proportions do not sum to one in modeled age ranges.")
  }

  if (any(!is.finite(unlist(sim[c("obs_catch", "ob_bts", "ob_ats", "obs_cpue", "ob_avo")])))) {
    stop("Simulated index/catch values include non-finite values.")
  }
  if (any(unlist(sim[c("obs_catch", "ob_bts", "ob_ats", "obs_cpue", "ob_avo")]) <= 0)) {
    stop("Simulated index/catch values must be positive.")
  }

  if (!is.null(report)) {
    report_checks <- c(
      length(report$pred_catch) == length(sim$obs_catch),
      length(report$eb_bts) == length(sim$ob_bts),
      length(report$eb_ats) == length(sim$ob_ats),
      length(report$pred_cpue) == length(sim$obs_cpue),
      length(report$pred_avo) == length(sim$ob_avo),
      all(dim(report$phat_fsh) == dim(sim$oac_fsh)),
      all(dim(report$phat_bts) == dim(sim$oac_bts)),
      all(dim(report$phat_ats) == dim(sim$oac_ats))
    )
    if (!all(report_checks)) {
      stop("Simulated data dimensions do not align with base-report fitted quantities.")
    }
  }

  TRUE
}

simulate_rtmb_datasets <- function(report, data, nsim = 1L, seed = NULL) {
  nsim <- as.integer(nsim)
  if (nsim < 1L) stop("nsim must be at least 1.")
  seeds <- if (is.null(seed)) rep(list(NULL), nsim) else as.list(as.integer(seed) + seq_len(nsim) - 1L)
  lapply(seeds, function(s) {
    sim <- prepare_simulated_estimation_data(simulate_rtmb_dataset(report, data, seed = s))
    validate_simulated_estimation_data(sim, report = report)
    sim
  })
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1] %||% "")
is_main_script <- nzchar(file_arg) && identical(basename(file_arg), "simulate_base_datasets.R")

if (is_main_script && !interactive()) {
  rtmb_root <- normalizePath(getwd(), mustWork = TRUE)
  base_path <- file.path(rtmb_root, "analysis", "output", "base.rds")
  out_path <- Sys.getenv(
    "OUTPUT",
    file.path(rtmb_root, "analysis", "output", "simulated", "base_simulated_datasets.rds")
  )
  nsim <- as.integer(Sys.getenv("NSIM", "10"))
  seed <- as.integer(Sys.getenv("SEED", "123"))

  pollock_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
  if (is.na(pollock_root) || !nzchar(pollock_root)) {
    pollock_root <- Sys.getenv("POLLOCK_BASE", unset = dirname(rtmb_root))
  }
  pollock_root <- normalizePath(pollock_root, mustWork = TRUE)

  base::source(file.path(rtmb_root, "R", "utilities.R"))
  data <- Get_Data()

  saved <- readRDS(base_path)
  sims <- simulate_rtmb_datasets(saved$report, data, nsim = nsim, seed = seed)
  validation <- vapply(sims, validate_simulated_estimation_data, logical(1), report = saved$report)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      simulations = sims,
      validation = validation,
      metadata = list(
        nsim = nsim,
        seed = seed,
        base_path = base_path,
        estimation_ready = all(validation),
        created = Sys.time()
      )
    ),
    out_path
  )
  cat("Wrote simulated data sets:", out_path, "\n")
}
