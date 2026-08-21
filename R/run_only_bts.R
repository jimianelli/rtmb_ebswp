# Fit the RTMB-ADMB bridge with BTS retained and ATS/AVO survey indices excluded.
# The ATS age-composition likelihood remains active; only ATS biomass, ATS age-1,
# and AVO index likelihood components are removed from the objective.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  config_path <- file.path(script_dir, "config.R")
} else {
  config_path <- file.path("R", "config.R")
}

source(config_path)

# config.R clears its sourcing environment as part of the historical bridge
# setup, so resolve the repository path again for output handling.
rtmb_dir <- normalizePath(getwd(), mustWork = TRUE)
admb_rep_path <- normalizePath(file.path(admb_run_dir, "pm.rep"), mustWork = TRUE)
admb_par_path <- normalizePath(file.path(admb_run_dir, "pm.par"), mustWork = TRUE)

fixed_params <- names(map_obj)[vapply(
  map_obj,
  function(x) all(is.na(x)),
  logical(1)
)]

parms$include_ats_index <- 0
parms$include_ats_age1_index <- 0
parms$include_avo_index <- 0
map_obj <- create_map_from_par(
  parms, parms,
  exact_names = c(fixed_params, "include_ats_index", "include_ats_age1_index", "include_avo_index"),
  exclude_patterns = "xxx"
)
data$return_nll_only <- 1

obj_only_bts <- MakeADFun(rpm, parms, map = map_obj)

start_par <- obj_only_bts$par
default_output_path <- if (identical(bridge_case, "corrected_full_age_bts")) {
  file.path(rtmb_dir, "analysis", "output", "corrected_full_age_bts", "only_bts.rds")
} else {
  file.path(rtmb_dir, "analysis", "output", "only_bts.rds")
}
output_path <- Sys.getenv("RTMB_ONLY_BTS_OUTPUT", unset = default_output_path)
if (file.exists(output_path)) {
  previous <- readRDS(output_path)
  same_configuration <- identical(
    previous$metadata$configuration,
    bridge_case
  )
  if (same_configuration && !is.null(previous$fit$par) &&
      length(previous$fit$par) == length(start_par)) {
    start_par <- previous$fit$par
  }
}

fit_only_bts <- nlminb(
  start_par,
  obj_only_bts$fn,
  obj_only_bts$gr,
  control = list(iter.max = 3000, eval.max = 3000)
)

first_gradient <- max(abs(obj_only_bts$gr(fit_only_bts$par)), na.rm = TRUE)
optimization_passes <- 1L
if (is.finite(first_gradient) && first_gradient > 1e-3) {
  fit_restart <- nlminb(
    fit_only_bts$par,
    obj_only_bts$fn,
    obj_only_bts$gr,
    control = list(
      iter.max = 3000,
      eval.max = 5000,
      rel.tol = 1e-12,
      x.tol = 1e-10
    )
  )
  if (is.finite(fit_restart$objective) &&
      fit_restart$objective <= fit_only_bts$objective + 1e-8) {
    fit_only_bts <- fit_restart
    optimization_passes <- 2L
  }
}

obj_only_bts$fn(fit_only_bts$par)
max_gradient <- max(abs(obj_only_bts$gr(fit_only_bts$par)), na.rm = TRUE)
fitted_parms <- obj_only_bts$env$parList(fit_only_bts$par)
data$return_nll_only <- 0
rtmb_report <- rpm(fitted_parms)

output_dir <- dirname(output_path)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
base_file <- Sys.getenv(
  "RTMB_ONLY_BTS_BASE_FILE",
  unset = if (identical(bridge_case, "corrected_full_age_bts")) {
    file.path(rtmb_dir, "analysis", "output", bridge_case, "rtmb_base.rds")
  } else {
    file.path(rtmb_dir, "analysis", "output", "base.rds")
  }
)

saveRDS(
  list(
    report = rtmb_report,
    fit = fit_only_bts,
    metadata = list(
      model = "Only BTS",
      configuration = bridge_case,
      bts_comp_normalization = bts_comp_normalization,
      base_file = base_file,
      base_md5 = if (file.exists(base_file)) {
        unname(tools::md5sum(base_file))
      } else {
        NA_character_
      },
      created = Sys.time(),
      admb_rep = admb_rep_path,
      admb_par = admb_par_path,
      max_gradient = max_gradient,
      optimization_passes = optimization_passes,
      excluded_likelihoods = c("ATS index", "ATS age-1 index", "AVO index")
    )
  ),
  output_path
)

cat("Wrote Only BTS RTMB output:", output_path, "\n")
cat("Convergence:", fit_only_bts$convergence, "\n")
cat("Objective:", fit_only_bts$objective, "\n")
cat("Maximum gradient:", max_gradient, "\n")
