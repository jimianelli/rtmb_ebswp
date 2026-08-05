#!/usr/bin/env Rscript

# Fit the 2D AR1 fishery selectivity form with the default AR1 penalty weight.
#
# This is a comparison companion to R/run_fishery_selectivity_forms.R, whose
# Form 5 scenario deliberately relaxes sel_tv_ar1_weight_fsh to 0.25. This
# script keeps sel_tv_ar1_weight_fsh = 1.0 and saves a separate result.

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

suppressPackageStartupMessages({
  library(RTMB)
  library(ebswp)
})

rtmb_root <- get_rtmb_root()
pollock_root <- get_pollock_root(rtmb_root)
out_dir <- file.path(rtmb_root, "analysis", "output", "fishery_sel_forms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

form_id <- 5L
data <- model_env$Get_Data()
data$sam_bts <- floor(pm$sam_bts)
data$fishery_sel_form <- form_id
data$return_nll_only <- 1
data$sel_tv_ar1_weight_fsh <- 1.0

parms <- model_env$read_pars(admb_par_path)
parms$steepness <- 0.67
parms <- model_env$add_fishery_selectivity_parameters(parms, data)

set.seed(123)
nyrs <- as.integer(data$endyr - data$styr + 1L)
nages <- as.integer(data$nages)
parms$log_sel_tv_ar1_sigma_fsh <- log(1.2)
parms$sel_tv_ar1_rho_fsh <- c(0, 0)
parms$sel_tv_ar1_fsh <- matrix(rnorm(nyrs * nages, sd = 0.05), nrow = nyrs, ncol = nages)

start_file <- Sys.getenv(
  "START_FILE",
  file.path(out_dir, "fishery_sel_form_5.rds")
)
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

cat("\n--- Fitting fishery_sel_form=5 with sel_tv_ar1_weight_fsh=1.0 ---\n")
t0 <- Sys.time()
fit <- nlminb(
  obj$par, obj$fn, obj$gr,
  control = list(eval.max = 5000, iter.max = 3000)
)
t1 <- Sys.time()

fitted_parms <- obj$env$parList(fit$par)
model_env$data$return_nll_only <- 0
report <- rpm(fitted_parms)
if (is.list(report) && !is.null(report$rtmb)) {
  report <- report$rtmb
}

res <- list(
  form = form_id,
  scenario = "default_weight",
  label = "2D AR1 year-age (default penalty weight)",
  objective = fit$objective,
  evaluated_total = report$tot_like %||% NA_real_,
  convergence = fit$convergence,
  message = fit$message,
  max_gradient = max(abs(obj$gr(fit$par)), na.rm = TRUE),
  seconds = as.numeric(difftime(t1, t0, units = "secs")),
  sel_tv_ar1_weight_fsh = data$sel_tv_ar1_weight_fsh,
  start_file = start_file,
  fixed_parameters = fitted_parms,
  report = report
)

out_file <- file.path(out_dir, "fishery_sel_form_5_default_weight.rds")
saveRDS(res, out_file)

summary <- data.frame(
  fishery_sel_form = form_id,
  scenario = res$scenario,
  label = res$label,
  objective = res$objective,
  evaluated_total = res$evaluated_total,
  convergence = res$convergence,
  max_gradient = res$max_gradient,
  seconds = res$seconds,
  sel_tv_ar1_weight_fsh = res$sel_tv_ar1_weight_fsh,
  start_file = res$start_file
)
utils::write.csv(
  summary,
  file.path(out_dir, "summary_form_5_default_weight.csv"),
  row.names = FALSE,
  na = ""
)

print(summary)
cat("\nWrote result to:", out_file, "\n")
