#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

script_path <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_path, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "scripts/check_fishery_selectivity_forms.R"
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo, "R", "model_funs.R"))

make_fishery_sel_spline_basis <- function(nages, nbasis = 6L, degree = 3L) {
  ages <- 0.5 + seq_len(nages)
  as.matrix(splines::bs(ages, df = nbasis, degree = degree, intercept = TRUE))
}

styr <- 2000L
endyr <- 2006L
nages <- 15L
nsel <- 12L
nyrs <- endyr - styr + 1L
yrs_ch_fsh <- styr:(endyr - 1L)

basis <- make_fishery_sel_spline_basis(nages, nbasis = 6)
common <- list(
  nsel = nsel,
  stsel = styr,
  endyr = endyr,
  nages = nages,
  coffs = seq(-2, 0.2, length.out = nsel),
  sel_devs = matrix(0.01, nrow = length(yrs_ch_fsh), ncol = nsel),
  yrs_ch_fsh = yrs_ch_fsh,
  sel_logistic_fsh = c(log(5), log(1.5)),
  sel_double_logistic_fsh = log(c(1.5, 3, 2.5)),
  sel_richards_fsh = c(log(4), log(1), log(1), log(5), log(0.75), log(1)),
  sel_spline_fsh = seq(-0.5, 0.4, length.out = ncol(basis)),
  fishery_sel_spline_basis = basis,
  sel_tv_ar1_fsh = matrix(rnorm(nyrs * nages, sd = 0.05), nrow = nyrs, ncol = nages)
)

for (form in 0:5) {
  out <- do.call(compute_selectivity_fsh_forms, c(list(fishery_sel_form = form), common))
  stopifnot(is.matrix(out$log_sel), all(dim(out$log_sel) == c(nyrs, nages)))
  stopifnot(all(is.finite(out$log_sel)), all(is.finite(exp(out$log_sel))))
  centered <- apply(exp(out$log_sel), 1, mean)
  stopifnot(max(abs(centered - 1)) < 1e-8)
  old_age_range <- 11:nages
  old_age_diff <- sweep(out$log_sel[, old_age_range, drop = FALSE], 1, out$log_sel[, 11], "-")
  stopifnot(max(abs(old_age_diff)) < 1e-8)
}

message("Fishery selectivity form generator check passed for forms 0:5 with ages 11+ tied.")
