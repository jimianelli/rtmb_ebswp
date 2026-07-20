# Run likelihood profiles for the default EBS pollock RTMB bridge model.
#
# Environment variables:
#   PROFILE_PARAMETERS  comma-separated names; repeated names use name[n]
#   PROFILE_POINTS      odd number of grid points (default 17 for log_avgrec)
#   PROFILE_HALF_WIDTH  half-width on fitted scale (default 2 for log_avgrec)
#   PROFILE_MODE        reopt or slice (default reopt)
#   PROFILE_MAX_EVAL    optimizer evaluation limit (default 5000)
#   PROFILE_GRADIENT_TOL maximum acceptable base gradient (default 0.002)
#   PROFILE_ALLOW_NONCONVERGED true to override the base-fit quality check

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath(file.path("R", "run_likelihood_profiles.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
config_path <- file.path(repo_root, "R", "config.R")

options(MakeADFun.silent = TRUE)
source(config_path)
repo_root <- rtmb_dir
source(file.path(repo_root, "R", "profile_components.R"))

parameters <- trimws(strsplit(
  Sys.getenv("PROFILE_PARAMETERS", "log_avgrec"), ",", fixed = TRUE
)[[1]])
default_avgrec_profile <- identical(parameters, "log_avgrec")
n_points <- as.integer(Sys.getenv(
  "PROFILE_POINTS", if (default_avgrec_profile) "17" else "15"
))
half_width <- as.numeric(Sys.getenv(
  "PROFILE_HALF_WIDTH", if (default_avgrec_profile) "2" else "0.35"
))
mode <- Sys.getenv("PROFILE_MODE", "reopt")
max_eval <- as.integer(Sys.getenv("PROFILE_MAX_EVAL", "5000"))
gradient_tolerance <- as.numeric(Sys.getenv("PROFILE_GRADIENT_TOL", "0.002"))
allow_nonconverged <- tolower(Sys.getenv("PROFILE_ALLOW_NONCONVERGED", "false")) %in%
  c("1", "true", "yes")

if (length(parameters) < 1L || any(!nzchar(parameters))) {
  stop("PROFILE_PARAMETERS must contain at least one parameter name.", call. = FALSE)
}
if (is.na(n_points) || n_points < 3L || n_points %% 2L == 0L) {
  stop("PROFILE_POINTS must be an odd integer of at least 3.", call. = FALSE)
}
if (!is.finite(half_width) || half_width <= 0) {
  stop("PROFILE_HALF_WIDTH must be a positive number.", call. = FALSE)
}
if (!mode %in% c("reopt", "slice")) {
  stop("PROFILE_MODE must be either `reopt` or `slice`.", call. = FALSE)
}
if (is.na(max_eval) || max_eval < 1L) {
  stop("PROFILE_MAX_EVAL must be a positive integer.", call. = FALSE)
}
if (!is.finite(gradient_tolerance) || gradient_tolerance <= 0) {
  stop("PROFILE_GRADIENT_TOL must be a positive number.", call. = FALSE)
}

parameter_base_names <- sub("\\[[0-9]+\\]$", "", parameters)
missing_from_model <- setdiff(parameter_base_names, unique(names(obj$par)))
if (length(missing_from_model) > 0L) {
  unavailable <- setdiff(missing_from_model, names(parms))
  if (length(unavailable) > 0L) {
    stop(
      "Profile parameter(s) not found in the RTMB parameter list: ",
      paste(unavailable, collapse = ", "),
      call. = FALSE
    )
  }

  fixed_params <- setdiff(fixed_params, missing_from_model)
  map_obj <- create_map_from_par(
    parms, parms,
    exact_names = fixed_params,
    exclude_patterns = "xxx"
  )
  obj <- MakeADFun(rpm, parms, map = map_obj)
  message(
    "Released fixed parameter(s) for profiling: ",
    paste(missing_from_model, collapse = ", ")
  )
}

message("Fitting the base RTMB model before profiling...")
fit <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = max_eval, eval.max = max_eval)
)
base_gradient <- max(abs(obj$gr(fit$par)))
base_polish_steps <- 0L
while (is.finite(base_gradient) &&
       base_gradient > gradient_tolerance &&
       base_polish_steps < 3L) {
  previous_gradient <- base_gradient
  polished_fit <- nlminb(
    start = fit$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = max_eval, eval.max = max_eval)
  )
  polished_gradient <- max(abs(obj$gr(polished_fit$par)))
  if (polished_fit$objective <= fit$objective + 1e-8) {
    fit <- polished_fit
    base_gradient <- polished_gradient
  }
  base_polish_steps <- base_polish_steps + 1L
  if (base_gradient >= previous_gradient) {
    break
  }
}
message("Base objective: ", signif(fit$objective, 10))
message("Base maximum absolute gradient: ", signif(base_gradient, 6))
if (fit$convergence != 0L || base_gradient > gradient_tolerance) {
  quality_message <- paste0(
    "Base fit did not meet the profile quality check: convergence code ",
    fit$convergence, "; maximum absolute gradient ", signif(base_gradient, 6),
    " (limit ", gradient_tolerance, ")."
  )
  if (!allow_nonconverged) {
    stop(
      quality_message,
      " Increase PROFILE_MAX_EVAL or set PROFILE_ALLOW_NONCONVERGED=true for diagnostic use.",
      call. = FALSE
    )
  }
  warning(quality_message)
}

profile_fun <- if (mode == "reopt") {
  profile_components_reopt
} else {
  profile_components_slice
}

results <- lapply(parameters, function(parameter) {
  target <- profile_resolve_parameter(fit$par, parameter)
  center <- fit$par[target$index]
  grid <- seq(center - half_width, center + half_width, length.out = n_points)
  profile_args <- list(
    obj = obj, name = parameter, values = grid, par = fit$par
  )
  if (mode == "reopt") {
    profile_args$control <- list(iter.max = max_eval, eval.max = max_eval)
    profile_args$gradient_tolerance <- gradient_tolerance
  }
  do.call(profile_fun, profile_args)
})
names(results) <- parameters

output_dir <- file.path(repo_root, "analysis", "output", "profiles")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- list(
  created = Sys.time(),
  mode = mode,
  parameters = parameters,
  points = n_points,
  half_width = half_width,
  gradient_tolerance = gradient_tolerance,
  base_fit = fit,
  base_max_gradient = base_gradient,
  base_polish_steps = base_polish_steps
)
saveRDS(
  list(metadata = metadata, profiles = results),
  file.path(output_dir, "likelihood_profiles.rds")
)

for (parameter in parameters) {
  file_label <- gsub("[^A-Za-z0-9_-]+", "_", parameter)
  write.csv(
    results[[parameter]],
    file.path(output_dir, paste0("profile_", file_label, ".csv")),
    row.names = FALSE
  )
  plot <- plot_profile_components(results[[parameter]])
  ggplot2::ggsave(
    file.path(output_dir, paste0("profile_", file_label, ".png")),
    plot = plot,
    width = 10,
    height = 8,
    dpi = 150
  )
}

message("Wrote profile outputs to: ", output_dir)
