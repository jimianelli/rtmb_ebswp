# Run strict RTMBprof component profiles for candidate abundance-scale
# parameters in the default EBS pollock RTMB bridge model.
#
# Environment variables:
#   RTMBPROF_PARAMETERS comma-separated parameter names
#   RTMBPROF_YTOL       target objective rise (default 2.2)
#   RTMBPROF_YSTEP      adaptive objective resolution (default 0.2)
#   RTMBPROF_H          initial focal-parameter step (default 0.02)
#   RTMBPROF_MAXIT      maximum points in each direction (default 35)
#   RTMBPROF_MAX_EVAL   nuisance optimizer limit (default 5000)
#   RTMBPROF_REL_TOL    nuisance optimizer relative tolerance (default 1e-8)

if (!requireNamespace("RTMBprof", quietly = TRUE)) {
  stop("Install RTMBprof before running this script.", call. = FALSE)
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath(file.path("R", "run_rtmbprof_stage1.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

parameters <- trimws(strsplit(
  Sys.getenv(
    "RTMBPROF_PARAMETERS",
    "log_q_ats"
  ),
  ",",
  fixed = TRUE
)[[1]])
ytol <- as.numeric(Sys.getenv("RTMBPROF_YTOL", "2.2"))
ystep <- as.numeric(Sys.getenv("RTMBPROF_YSTEP", "0.2"))
h <- as.numeric(Sys.getenv("RTMBPROF_H", "0.02"))
maxit <- as.integer(Sys.getenv("RTMBPROF_MAXIT", "35"))
max_eval <- as.integer(Sys.getenv("RTMBPROF_MAX_EVAL", "5000"))
rel_tol <- as.numeric(Sys.getenv("RTMBPROF_REL_TOL", "1e-8"))

if (length(parameters) < 1L || any(!nzchar(parameters))) {
  stop("RTMBPROF_PARAMETERS must name at least one parameter.", call. = FALSE)
}
if (!is.finite(ytol) || ytol <= 0 ||
    !is.finite(ystep) || ystep <= 0 ||
    !is.finite(h) || h <= 0 ||
    is.na(maxit) || maxit < 1L ||
    is.na(max_eval) || max_eval < 1L ||
    !is.finite(rel_tol) || rel_tol <= 0) {
  stop("Profile controls must be finite and positive.", call. = FALSE)
}

options(MakeADFun.silent = TRUE)
rtmb_env <- new.env(parent = globalenv())
rtmb_env$rm <- function(...) invisible(NULL)
rtmb_env$source <- function(file, ...) {
  base::source(file, local = parent.frame(), ...)
}
base::source(file.path(repo_root, "R", "config.R"), local = rtmb_env)
base::source(file.path(repo_root, "R", "profile_components.R"))

obj <- rtmb_env$obj
obj$env$tracemgc <- FALSE
missing_parameters <- setdiff(parameters, unique(names(obj$par)))
if (length(missing_parameters) > 0L) {
  stop(
    "Candidate parameter(s) are not active in the default model: ",
    paste(missing_parameters, collapse = ", "),
    call. = FALSE
  )
}

fit <- stats::nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = max_eval, eval.max = max_eval)
)
for (polish in seq_len(3L)) {
  max_gradient <- max(abs(obj$gr(fit$par)))
  if (fit$convergence == 0L && max_gradient <= 0.002) {
    break
  }
  candidate <- stats::nlminb(
    start = fit$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = max_eval, eval.max = max_eval)
  )
  if (candidate$objective <= fit$objective + 1e-8) {
    fit <- candidate
  }
}
invisible(obj$fn(fit$par))
base_max_gradient <- max(abs(obj$gr(fit$par)))
if (fit$convergence != 0L || base_max_gradient > 0.002) {
  stop(
    "Base fit failed the Stage 1 quality check: convergence code ",
    fit$convergence,
    "; maximum absolute gradient ",
    signif(base_max_gradient, 6),
    ".",
    call. = FALSE
  )
}

output_dir <- file.path(repo_root, "analysis", "output", "rtmbprof_stage1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

profiles <- vector("list", length(parameters))
names(profiles) <- parameters
for (parameter in parameters) {
  message("Profiling candidate scale parameter: ", parameter)
  profiles[[parameter]] <- RTMBprof::profile_components(
    obj = obj,
    name = parameter,
    components = pollock_profile_components,
    objective = "joint",
    h = h,
    ytol = ytol,
    ystep = ystep,
    maxit = maxit,
    trace = FALSE,
    control = list(
      step.min = 0.001,
      iter.max = max_eval,
      eval.max = max_eval,
      rel.tol = rel_tol
    )
  )
  saveRDS(
    profiles[[parameter]],
    file.path(
      output_dir,
      paste0("checkpoint_", gsub("[^A-Za-z0-9_-]+", "_", parameter), ".rds")
    )
  )
}

scale_diagnostics <- lapply(profiles, function(profile) {
  constrained_parameters <- attr(profile, "parameters")
  reports <- lapply(constrained_parameters, obj$report)
  data.frame(
    parameter = attr(profile, "target"),
    focal_value = profile[[attr(profile, "target")]],
    objective = profile$value,
    terminal_SSB = vapply(
      reports,
      function(report) tail(report$SSB, 1L),
      numeric(1)
    ),
    mean_SSB = vapply(
      reports,
      function(report) mean(report$SSB),
      numeric(1)
    ),
    Bzero = vapply(reports, `[[`, numeric(1), "Bzero"),
    q_bts = vapply(reports, `[[`, numeric(1), "q_bts"),
    closure_error = profile$closure_error,
    convergence = profile$convergence,
    max_nuisance_gradient = profile$max_nuisance_gradient,
    check.names = FALSE
  )
})

profile_summaries <- lapply(parameters, function(parameter) {
  profile <- profiles[[parameter]]
  diagnostics <- scale_diagnostics[[parameter]]
  focal <- diagnostics$focal_value
  delta_objective <- diagnostics$objective - min(diagnostics$objective)
  fitted_value <- focal[which.min(abs(delta_objective))]
  left <- focal < fitted_value
  right <- focal > fitted_value

  data.frame(
    parameter = parameter,
    points = nrow(diagnostics),
    fitted_value = fitted_value,
    focal_min = min(focal),
    focal_max = max(focal),
    max_delta_objective = max(delta_objective),
    threshold_crossed_left = any(delta_objective[left] >= 1.92),
    threshold_crossed_right = any(delta_objective[right] >= 1.92),
    max_abs_closure_error = max(abs(diagnostics$closure_error)),
    max_nuisance_gradient = max(diagnostics$max_nuisance_gradient),
    nonzero_convergence_codes = sum(diagnostics$convergence != 0L),
    terminal_SSB_min = min(diagnostics$terminal_SSB),
    terminal_SSB_max = max(diagnostics$terminal_SSB),
    mean_SSB_min = min(diagnostics$mean_SSB),
    mean_SSB_max = max(diagnostics$mean_SSB),
    Bzero_min = min(diagnostics$Bzero),
    Bzero_max = max(diagnostics$Bzero),
    q_bts_min = min(diagnostics$q_bts),
    q_bts_max = max(diagnostics$q_bts),
    check.names = FALSE
  )
})
profile_summary <- do.call(rbind, profile_summaries)

invisible(obj$fn(fit$par))
output_file <- file.path(output_dir, "candidate_scale_profiles.rds")
saveRDS(
  list(
    created = Sys.time(),
    parameters = parameters,
    base_fit = fit,
    base_max_gradient = base_max_gradient,
    profiles = profiles,
    scale_diagnostics = scale_diagnostics,
    profile_summary = profile_summary
  ),
  output_file
)

for (parameter in parameters) {
  file_label <- gsub("[^A-Za-z0-9_-]+", "_", parameter)
  utils::write.csv(
    scale_diagnostics[[parameter]],
    file.path(output_dir, paste0(file_label, "_scale_diagnostics.csv")),
    row.names = FALSE
  )
}
utils::write.csv(
  profile_summary,
  file.path(output_dir, "profile_summary.csv"),
  row.names = FALSE
)

message("Wrote strict RTMBprof Stage 1 output to: ", output_file)
