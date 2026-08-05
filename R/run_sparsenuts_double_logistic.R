#!/usr/bin/env Rscript

# Run SparseNUTS for the time-varying double-logistic fishery selectivity form.
#
# Usage:
#   Rscript R/run_sparsenuts_double_logistic.R
#
# Optional environment variables:
#   SPARSENUTS_CHAINS=6
#   SPARSENUTS_CORES=6
#   SPARSENUTS_ITER=1000
#   SPARSENUTS_WARMUP=500
#   SPARSENUTS_SEED=123
#
# Outputs:
#   analysis/output/sparsenuts/fishery_sel_forms/rtmb_ebswp_sparsenuts_form_2.rds
#   analysis/output/sparsenuts/fishery_sel_forms/sparsenuts_diagnostics_form_2.rds
#   analysis/output/sparsenuts/fishery_sel_forms/sparsenuts_slow_parameters_form_2.csv
#   analysis/output/sparsenuts/fishery_sel_forms/figures/sparsenuts_pairs_tmbfit_form_2.png

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

get_rtmb_root <- function() {
  env_root <- Sys.getenv("RTMB_EBSWP_ROOT", unset = NA_character_)
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

get_pollock_root <- function(rtmb_root) {
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
  candidates <- c(rtmb_root, dirname(rtmb_root), file.path(dirname(rtmb_root), "pollock"))
  for (cand in candidates) {
    if (has_bridge(cand)) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  stop("Cannot locate pollock bridge inputs. Set POLLOCK_ROOT or use bundled admb/runs files.")
}

write_csv_base <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

select_slow_parameters <- function(fit, n = 6L) {
  mon <- as.data.frame(fit$monitor)
  if (!nrow(mon)) {
    return(mon)
  }
  if (!"variable" %in% names(mon)) {
    mon$variable <- rownames(mon)
  }
  if (!"rhat" %in% names(mon)) mon$rhat <- NA_real_
  if (!"ess_bulk" %in% names(mon)) mon$ess_bulk <- NA_real_

  rhat_rank <- ifelse(is.finite(mon$rhat), -mon$rhat, Inf)
  ess_rank <- ifelse(is.finite(mon$ess_bulk), mon$ess_bulk, Inf)
  mon[order(rhat_rank, ess_rank), , drop = FALSE][seq_len(min(n, nrow(mon))), , drop = FALSE]
}

save_png <- function(file, expr, width = 1800, height = 1600, res = 180) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = width, height = height, res = res)
  ok <- tryCatch(
    {
      force(expr)
      TRUE
    },
    error = function(e) {
      message("Plot failed for ", basename(file), ": ", conditionMessage(e))
      FALSE
    }
  )
  grDevices::dev.off()
  if (!ok) {
    unlink(file)
  }
  invisible(ok)
}

suppressPackageStartupMessages({
  library(RTMB)
  library(ebswp)
  library(SparseNUTS)
})

rtmb_root <- get_rtmb_root()
pollock_root <- get_pollock_root(rtmb_root)

out_dir <- file.path(rtmb_root, "analysis", "output", "sparsenuts", "fishery_sel_forms")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

chains <- as.integer(Sys.getenv("SPARSENUTS_CHAINS", "6"))
cores <- as.integer(Sys.getenv("SPARSENUTS_CORES", "6"))
iter <- as.integer(Sys.getenv("SPARSENUTS_ITER", "1000"))
warmup <- as.integer(Sys.getenv("SPARSENUTS_WARMUP", as.character(floor(iter / 2))))
num_samples <- max(1L, iter - warmup)
seed <- as.integer(Sys.getenv("SPARSENUTS_SEED", "123"))
form_id <- 2L

model_env <- new.env(parent = globalenv())
source(file.path(rtmb_root, "R", "utilities.R"), local = model_env)
source(file.path(rtmb_root, "R", "model_funs.R"), local = model_env)
source(file.path(rtmb_root, "R", "Rpm.R"), local = model_env)

admb_rep_path <- file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.rep")
if (!file.exists(admb_rep_path)) {
  admb_rep_path <- file.path(pollock_root, "runs", "for_rtmb", "pm.rep")
}
if (!file.exists(admb_rep_path)) {
  stop("ADMB rep file not found.")
}

pm <- read_rep(admb_rep_path)
pm$phat_ats <- pm$phat_ats[, 2:16]
pm$phat_bts <- pm$phat_bts[, 2:16]
pm$sel_ats <- pm$sel_ats[31:61, ]
pm$sel_bts <- pm$sel_bts[19:61, ]
pm$phat_fsh <- pm$phat_fsh[, 2:16]
pm$bts_like <- pm$surv_like[1]
pm$ats_like <- pm$surv_like[2]
pm$ats_age1_like <- pm$surv_like[3]
pm$SSB <- pm$SSB[, 2]
model_env$pm <- pm

admb_par_path <- file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.par")
if (!file.exists(admb_par_path)) {
  admb_par_path <- file.path(pollock_root, "runs", "rtmb", "pm.par")
}
if (!file.exists(admb_par_path)) {
  stop("ADMB par file not found.")
}

data <- model_env$Get_Data()
data$sam_bts <- floor(pm$sam_bts)
data$fishery_sel_form <- form_id
data$return_nll_only <- 1

parms <- model_env$read_pars(admb_par_path)
parms$steepness <- 0.67
parms <- model_env$add_fishery_selectivity_parameters(parms, data)

start_file <- file.path(rtmb_root, "analysis", "output", "fishery_sel_forms", "fishery_sel_form_2.rds")
if (file.exists(start_file)) {
  start_saved <- readRDS(start_file)
  if (!is.null(start_saved$fixed_parameters)) {
    parms[names(start_saved$fixed_parameters)] <- start_saved$fixed_parameters
  }
}

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
  model_env$inactive_fishery_selectivity_parameters(form_id),
  "resid_temp_x1", "resid_temp_x2",
  "log_q_std_area", "bt_slope", "sigr", "steepness",
  "sel_coffs_bts",
  "sel_tv_ar1_weight_fsh",
  "M_dev",
  "log_a_II", "log_b_II", "log_a_II_vec", "log_b_II_vec",
  "log_rho", "log_resid_M",
  "log_alpha"
)

map_obj <- model_env$create_map_from_par(
  parms, parms,
  exact_names = fixed_params,
  exclude_patterns = "xxx"
)

rpm <- model_env$rpm
environment(rpm) <- model_env
model_env$data <- data
obj <- RTMB::MakeADFun(rpm, parms, map = map_obj)

lower <- rep(-Inf, length(obj$par))
upper <- rep(Inf, length(obj$par))
idx <- which(names(obj$par) == "sel_double_logistic_fsh")
if (length(idx) > 0L) {
  nyrs <- length(idx) / 3L
  lower[idx] <- rep(log(c(0.25, 1.0, 1.5)), each = nyrs)
  upper[idx] <- rep(log(c(5.0, 8.0, 25.0)), each = nyrs)
}

cat(
  "Running SparseNUTS for fishery_sel_form=2 with",
  chains, "chains,", cores, "cores,", iter, "iterations,",
  warmup, "warmup iterations.\n"
)

t0 <- Sys.time()
fit <- SparseNUTS::sample_snuts(
  obj,
  chains = chains,
  cores = cores,
  num_samples = num_samples,
  num_warmup = warmup,
  seed = seed,
  init = "last.par.best",
  metric = "unit",
  lower = lower,
  upper = upper,
  globals = list(data = data),
  model_name = "rtmb_ebswp_double_logistic"
)
t1 <- Sys.time()

metadata <- list(
  fishery_sel_form = form_id,
  label = "Time-varying double logistic (annual p1/p2/p3)",
  chains = chains,
  cores = cores,
  iter = iter,
  warmup = warmup,
  num_samples = num_samples,
  seed = seed,
  metric = "unit",
  init = "last.par.best",
  lower_upper = "double-logistic annual parameters constrained to optimization bounds",
  start_file = start_file,
  seconds = as.numeric(difftime(t1, t0, units = "secs")),
  created = Sys.time()
)
attr(fit, "rtmb_ebswp_sparsenuts") <- metadata

fit_file <- file.path(out_dir, "rtmb_ebswp_sparsenuts_form_2.rds")
saveRDS(fit, fit_file)
saveRDS(metadata, file.path(out_dir, "sparsenuts_metadata_form_2.rds"))

diagnostics <- tryCatch(
  SparseNUTS::check_snuts_diagnostics(fit, print = FALSE),
  error = function(e) list(error = conditionMessage(e))
)
saveRDS(diagnostics, file.path(out_dir, "sparsenuts_diagnostics_form_2.rds"))
writeLines(
  capture.output(str(diagnostics, max.level = 2)),
  file.path(out_dir, "sparsenuts_diagnostics_form_2.txt")
)

slow <- select_slow_parameters(fit, n = 6L)
write_csv_base(slow, file.path(out_dir, "sparsenuts_slow_parameters_form_2.csv"))

available <- dimnames(fit$samples)[[3]]
slow_names <- slow$variable
slow_idx <- match(slow_names, available)
slow_idx <- slow_idx[is.finite(slow_idx)]
if (!length(slow_idx)) {
  slow_idx <- seq_len(min(6L, dim(fit$samples)[3]))
}

save_png(
  file.path(fig_dir, "sparsenuts_pairs_tmbfit_form_2.png"),
  pairs(
    fit,
    pars = slow_idx,
    order = "orig",
    diag = "hist",
    point.col = grDevices::adjustcolor("black", alpha.f = 0.25),
    point.pch = 16,
    plot = TRUE
  ),
  width = 1800,
  height = 1800
)

save_png(
  file.path(fig_dir, "sparsenuts_marginals_slow_form_2.png"),
  SparseNUTS::plot_marginals(fit, pars = slow_idx, mfrow = c(3, 2)),
  width = 1800,
  height = 1400
)

save_png(
  file.path(fig_dir, "sparsenuts_sampler_params_form_2.png"),
  SparseNUTS::plot_sampler_params(fit, plot = TRUE),
  width = 1600,
  height = 1200
)

cat("Wrote SparseNUTS fit:", fit_file, "\n")
