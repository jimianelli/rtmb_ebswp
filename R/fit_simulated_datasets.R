#!/usr/bin/env Rscript

# Refit the RTMB bridge model to each simulated pseudo-data set and save
# recruitment and spawning biomass summaries.
#
# Usage:
#   MAXSIM=30 Rscript R/fit_simulated_datasets.R
#   FISHERY_SEL_FORM=5 OUTPUT_DIR=analysis/output/simulated/refits_2dar1 \
#     MAXSIM=30 Rscript R/fit_simulated_datasets.R
#
# Outputs:
#   analysis/output/simulated/refits/simulated_refits.rds
#   analysis/output/simulated/refits/simulated_refit_diagnostics.csv
#   analysis/output/simulated/refits/simulated_refit_timeseries.csv

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
  env_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
  if (is.na(env_root) || !nzchar(env_root)) {
    env_root <- Sys.getenv("POLLOCK_BASE", unset = NA_character_)
  }
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }
  normalizePath(dirname(rtmb_root), mustWork = TRUE)
}

write_csv_base <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

as_series <- function(x, sim, quantity) {
  data.frame(
    simulation = sim,
    quantity = quantity,
    year = as.integer(names(x) %||% seq_along(x)),
    value = as.numeric(x)
  )
}

as_selectivity_series <- function(x, sim, ages = 3:8) {
  x <- as.matrix(x)
  cols <- match(as.character(ages), colnames(x))
  if (anyNA(cols)) {
    cols <- ages
  }
  out <- expand.grid(
    year = as.integer(rownames(x) %||% seq_len(nrow(x))),
    age = ages,
    KEEP.OUT.ATTRS = FALSE
  )
  out$simulation <- sim
  out$value <- as.numeric(x[, cols, drop = FALSE])
  out[, c("simulation", "year", "age", "value")]
}

rtmb_dir <- get_rtmb_root()
pollock_root <- get_pollock_root(rtmb_dir)

suppressPackageStartupMessages({
  library(RTMB)
  library(ebswp)
})

model_env <- new.env(parent = globalenv())
source(file.path(rtmb_dir, "R", "utilities.R"), local = model_env)
source(file.path(rtmb_dir, "R", "model_funs.R"), local = model_env)
source(file.path(rtmb_dir, "R", "Rpm.R"), local = model_env)
source(file.path(rtmb_dir, "R", "simulate_base_datasets.R"), local = model_env)
rpm <- model_env$rpm

sim_file <- Sys.getenv(
  "SIM_FILE",
  file.path(rtmb_dir, "analysis", "output", "simulated", "base_simulated_datasets.rds")
)
fishery_sel_form <- as.integer(Sys.getenv("FISHERY_SEL_FORM", "0"))
fishery_sel_label <- Sys.getenv(
  "FISHERY_SEL_LABEL",
  switch(as.character(fishery_sel_form),
         "0" = "base",
         "1" = "logistic",
         "2" = "double_logistic",
         "3" = "richards",
         "4" = "spline",
         "5" = "2dar1",
         paste0("form_", fishery_sel_form))
)
use_comparison_recipe <- as.logical(as.integer(Sys.getenv(
  "USE_FISHERY_SEL_COMPARISON_RECIPE",
  if (fishery_sel_form == 5L) "1" else "0"
)))
out_dir <- Sys.getenv(
  "OUTPUT_DIR",
  if (fishery_sel_form == 0L) {
    file.path(rtmb_dir, "analysis", "output", "simulated", "refits")
  } else {
    file.path(rtmb_dir, "analysis", "output", "simulated", paste0("refits_", fishery_sel_label))
  }
)
maxsim <- as.integer(Sys.getenv("MAXSIM", "30"))
eval_max <- as.integer(Sys.getenv("EVAL_MAX", "3000"))
iter_max <- as.integer(Sys.getenv("ITER_MAX", "3000"))
start_file <- Sys.getenv(
  "START_FILE",
  if (fishery_sel_form == 0L) {
    file.path(rtmb_dir, "analysis", "output", "fishery_sel_forms", "fishery_sel_form_0.rds")
  } else {
    ""
  }
)

if (!file.exists(sim_file)) {
  stop("Missing simulated datasets: ", sim_file)
}

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

base_data <- model_env$Get_Data()
if (!is.null(pm)) {
  base_data$sam_bts <- floor(pm$sam_bts)
}
base_data$fishery_sel_form <- fishery_sel_form

parms <- model_env$read_pars(admb_par_path)
parms$steepness <- 0.67
parms <- model_env$add_fishery_selectivity_parameters(parms, base_data)
if (file.exists(start_file)) {
  start_saved <- readRDS(start_file)
  if (!is.null(start_saved$fixed_parameters)) {
    parms[names(start_saved$fixed_parameters)] <- start_saved$fixed_parameters
  }
}
if (isTRUE(use_comparison_recipe) && fishery_sel_form == 5L) {
  set.seed(123)
  nyrs <- as.integer(base_data$endyr - base_data$styr + 1L)
  nages <- as.integer(base_data$nages)
  parms$log_sel_tv_ar1_sigma_fsh <- log(1.2)
  parms$sel_tv_ar1_rho_fsh <- c(0, 0)
  parms$sel_tv_ar1_fsh <- matrix(rnorm(nyrs * nages, sd = 0.05), nrow = nyrs, ncol = nages)
  base_data$sel_tv_ar1_weight_fsh <- 0.25
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
  model_env$inactive_fishery_selectivity_parameters(fishery_sel_form),
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

sim_saved <- readRDS(sim_file)
simulations <- sim_saved$simulations
if (length(simulations) == 0) {
  stop("Simulation file contains no simulations.")
}
maxsim <- min(maxsim, length(simulations))

fit_one_simulation <- function(sim_id) {
  message("Fitting simulated dataset ", sim_id, " of ", maxsim)
  model_env$data <- model_env$prepare_simulated_estimation_data(simulations[[sim_id]])
  model_env$data$fishery_sel_form <- fishery_sel_form
  if (isTRUE(use_comparison_recipe) && fishery_sel_form == 5L) {
    model_env$data$sel_tv_ar1_weight_fsh <- base_data$sel_tv_ar1_weight_fsh
  }
  model_env$validate_simulated_estimation_data(model_env$data)

  obj_sim <- RTMB::MakeADFun(rpm, parms, map = map_obj)
  t0 <- Sys.time()
  fit <- nlminb(
    obj_sim$par,
    obj_sim$fn,
    obj_sim$gr,
    control = list(eval.max = eval_max, iter.max = iter_max)
  )
  t1 <- Sys.time()

  obj_sim$fn(fit$par)
  fitted_parms <- obj_sim$env$parList(fit$par)
  model_env$data$return_nll_only <- 0
  report <- rpm(fitted_parms)
  if (is.list(report) && !is.null(report$rtmb)) {
    report <- report$rtmb
  }
  if (is.null(report$recruitment) && !is.null(report$N)) {
    report$recruitment <- report$N[, 1]
  }
  model_env$data$return_nll_only <- 1

  list(
    simulation = sim_id,
    fit = fit,
    convergence = fit$convergence,
    message = fit$message,
    objective = fit$objective,
    max_gradient = max(abs(obj_sim$gr(fit$par)), na.rm = TRUE),
    seconds = as.numeric(difftime(t1, t0, units = "secs")),
    SSB = report$SSB,
    recruitment = report$recruitment,
    sel_fsh = report$sel_fsh,
    tot_like = report$tot_like
  )
}

results <- lapply(seq_len(maxsim), function(i) {
  tryCatch(
    fit_one_simulation(i),
    error = function(e) {
      message("Simulated dataset ", i, " failed: ", conditionMessage(e))
      list(
        simulation = i,
        fit = NULL,
        convergence = NA_integer_,
        message = conditionMessage(e),
        objective = NA_real_,
        max_gradient = NA_real_,
        seconds = NA_real_,
        SSB = NULL,
        recruitment = NULL,
        sel_fsh = NULL,
        tot_like = NA_real_
      )
    }
  )
})

diagnostics <- do.call(
  rbind,
  lapply(results, function(x) {
    data.frame(
      simulation = x$simulation,
      fishery_sel_form = fishery_sel_form,
      fishery_sel_label = fishery_sel_label,
      convergence = x$convergence,
      message = x$message %||% NA_character_,
      objective = x$objective,
      max_gradient = x$max_gradient,
      seconds = x$seconds,
      tot_like = x$tot_like
    )
  })
)

series <- do.call(
  rbind,
  unlist(
    lapply(results, function(x) {
      if (is.null(x$SSB) || is.null(x$recruitment)) {
        return(list(NULL))
      }
      list(
        as_series(x$SSB, x$simulation, "SSB"),
        as_series(x$recruitment, x$simulation, "recruitment")
      )
    }),
    recursive = FALSE
  )
)

selectivity_series <- do.call(
  rbind,
  lapply(results, function(x) {
    if (is.null(x$sel_fsh)) {
      return(NULL)
    }
    as_selectivity_series(x$sel_fsh, x$simulation, ages = 3:8)
  })
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    diagnostics = diagnostics,
    series = series,
    selectivity_series = selectivity_series,
    results = results,
    metadata = list(
      sim_file = sim_file,
      maxsim = maxsim,
      eval_max = eval_max,
      iter_max = iter_max,
      fishery_sel_form = fishery_sel_form,
      fishery_sel_label = fishery_sel_label,
      use_comparison_recipe = use_comparison_recipe,
      start_file = start_file,
      created = Sys.time()
    )
  ),
  file.path(out_dir, "simulated_refits.rds")
)
write_csv_base(diagnostics, file.path(out_dir, "simulated_refit_diagnostics.csv"))
if (!is.null(series)) {
  write_csv_base(series, file.path(out_dir, "simulated_refit_timeseries.csv"))
}
if (!is.null(selectivity_series)) {
  write_csv_base(selectivity_series, file.path(out_dir, "simulated_refit_fishery_selectivity.csv"))
}

cat("Wrote simulated refits:", file.path(out_dir, "simulated_refits.rds"), "\n")
