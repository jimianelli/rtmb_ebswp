#!/usr/bin/env Rscript

# Focused smoke test for the conditional two-stage Form 2 implementation.
# The accepted bridge itself is covered by tests/test_sept_2025_bridge.R.

suppressPackageStartupMessages(library(RTMB))
source("R/utilities.R")
source("R/model_funs.R")

nyrs <- 4L
nages <- 15L
parameters <- list(
  mean = log(c(1.5, 3, 4)),
  deviations = matrix(0, nrow = nyrs, ncol = 3)
)

objective <- function(par) {
  effective <- par$deviations
  for (i in seq_len(nrow(effective))) {
    for (j in seq_len(ncol(effective))) {
      effective[i, j] <- par$mean[j] + par$deviations[i, j]
    }
  }
  selectivity <- compute_selectivity_fsh_double_logistic(
    stsel = 2000L,
    endyr = 2000L + nyrs - 1L,
    nages = nages,
    parameters = effective,
    old_age_cap = 0L
  )$log_sel
  penalty <- selectivity_like_fsh_double_logistic(
    log_sel_fsh = selectivity,
    selCFsh = 1,
    domFish = 1,
    mean_parameters = par$mean,
    annual_deviations = par$deviations,
    process_cv = 0.20,
    old_age_cap = 0L
  )
  penalty$total + 0.01 * sum(selectivity^2)
}

obj <- RTMB::MakeADFun(objective, parameters, random = "deviations")
stopifnot(is.finite(obj$fn(obj$par)))
stopifnot(all(is.finite(obj$gr(obj$par))))

plain <- compute_selectivity_fsh_double_logistic(
  stsel = 2000L,
  endyr = 2000L,
  nages = nages,
  parameters = matrix(log(c(1.5, 3, 4)), 1L, 3L),
  old_age_cap = 0L
)$log_sel
stopifnot(nrow(plain) == 1L, ncol(plain) == nages)
stopifnot(abs(mean(exp(plain[1, ])) - 1) < 1e-10)
stopifnot(length(unique(round(plain[1, 11:15], 10))) > 1L)

cat("PASS: conditional Form 2 selectivity tapes and leaves old ages uncapped.\n")
