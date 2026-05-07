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
if (file.exists(file.path(rtmb_dir, "analysis", "output", "only_bts.rds"))) {
  previous <- readRDS(file.path(rtmb_dir, "analysis", "output", "only_bts.rds"))
  if (!is.null(previous$fit$par) && length(previous$fit$par) == length(start_par)) {
    start_par <- previous$fit$par
  }
}

fit_only_bts <- nlminb(
  start_par,
  obj_only_bts$fn,
  obj_only_bts$gr,
  control = list(iter.max = 3000, eval.max = 3000)
)

obj_only_bts$fn(fit_only_bts$par)
max_gradient <- max(abs(obj_only_bts$gr(fit_only_bts$par)), na.rm = TRUE)
fitted_parms <- obj_only_bts$env$parList(fit_only_bts$par)
data$return_nll_only <- 0
rtmb_report <- rpm(fitted_parms)

output_dir <- file.path(rtmb_dir, "analysis", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_path <- file.path(output_dir, "only_bts.rds")
saveRDS(
  list(
    report = rtmb_report,
    fit = fit_only_bts,
    metadata = list(
      model = "Only BTS",
      created = Sys.time(),
      admb_rep = admb_rep_path,
      admb_par = admb_par_path,
      max_gradient = max_gradient,
      excluded_likelihoods = c("ATS index", "ATS age-1 index", "AVO index")
    )
  ),
  output_path
)

cat("Wrote Only BTS RTMB output:", output_path, "\n")
cat("Convergence:", fit_only_bts$convergence, "\n")
cat("Objective:", fit_only_bts$objective, "\n")
cat("Maximum gradient:", max_gradient, "\n")
