# Evaluate the strict EBS pollock objective at saved SparseNUTS draws and
# compare the resulting posterior-objective summaries with the Stage 1
# constrained profile for log_q_ats.
#
# Environment variables:
#   RTMBPROF_MCMC_FILE  saved SparseNUTS tmbfit file
#   RTMBPROF_MAX_DRAWS maximum retained draws to evaluate (0 means all)
#   RTMBPROF_BINS      number of conditional-summary bins (default 20)

if (!requireNamespace("RTMBprof", quietly = TRUE)) {
  stop("Install RTMBprof before running this script.", call. = FALSE)
}
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Install posterior before running this script.", call. = FALSE)
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath(file.path("R", "run_rtmbprof_stage2.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_dir <- file.path(repo_root, "analysis", "output", "rtmbprof_stage2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mcmc_file <- Sys.getenv(
  "RTMBPROF_MCMC_FILE",
  file.path(
    repo_root,
    "analysis", "output", "sparsenuts",
    "rtmb_ebswp_sparsenuts_default.rds"
  )
)
stage1_file <- file.path(
  repo_root,
  "analysis", "output", "rtmbprof_stage1",
  "candidate_scale_profiles.rds"
)
max_draws <- as.integer(Sys.getenv("RTMBPROF_MAX_DRAWS", "0"))
bins <- as.integer(Sys.getenv("RTMBPROF_BINS", "20"))
if (is.na(max_draws) || max_draws < 0L) {
  stop("RTMBPROF_MAX_DRAWS must be zero or a positive integer.", call. = FALSE)
}
if (is.na(bins) || bins < 2L) {
  stop("RTMBPROF_BINS must be an integer of at least two.", call. = FALSE)
}
if (!file.exists(mcmc_file)) {
  stop("Saved SparseNUTS fit not found: ", mcmc_file, call. = FALSE)
}
if (!file.exists(stage1_file)) {
  stop("Run R/run_rtmbprof_stage1.R before Stage 2.", call. = FALSE)
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

fit <- readRDS(mcmc_file)
required_fit_fields <- c("samples", "warmup", "monitor")
missing_fit_fields <- setdiff(required_fit_fields, names(fit))
if (length(missing_fit_fields) > 0L) {
  stop(
    "Saved SparseNUTS fit is missing: ",
    paste(missing_fit_fields, collapse = ", "),
    call. = FALSE
  )
}
samples <- fit$samples
if (!is.array(samples) || length(dim(samples)) != 3L) {
  stop("fit$samples must be iterations by chains by variables.", call. = FALSE)
}
warmup <- as.integer(fit$warmup)
if (length(warmup) != 1L || is.na(warmup) ||
    warmup < 0L || warmup >= dim(samples)[1]) {
  stop("The saved SparseNUTS warmup metadata is invalid.", call. = FALSE)
}
sample_names <- dimnames(samples)[[3]]
if (is.null(sample_names) || any(!nzchar(sample_names))) {
  stop("The saved SparseNUTS samples lack parameter names.", call. = FALSE)
}

# SparseNUTS retains warmup at the beginning of each chain. RTMBprof's generic
# array adapter cannot infer tmbfit warmup metadata, so remove it explicitly.
post_warmup <- samples[
  seq.int(warmup + 1L, dim(samples)[1]),
  ,
  ,
  drop = FALSE
]
parameter_names <- fit$par_names
if (!is.character(parameter_names) ||
    length(parameter_names) != length(obj$par) ||
    anyNA(parameter_names) || any(!nzchar(parameter_names)) ||
    anyDuplicated(parameter_names)) {
  stop(
    "fit$par_names must uniquely map every fitted RTMB parameter.",
    call. = FALSE
  )
}
missing_parameters <- setdiff(parameter_names, sample_names)
if (length(missing_parameters) > 0L) {
  stop(
    "Posterior draws are missing fitted parameter(s): ",
    paste(missing_parameters, collapse = ", "),
    call. = FALSE
  )
}
post_warmup <- post_warmup[, , parameter_names, drop = FALSE]

total_draws <- dim(post_warmup)[1] * dim(post_warmup)[2]
if (max_draws > 0L && max_draws < total_draws) {
  draws_per_chain <- max(1L, floor(max_draws / dim(post_warmup)[2]))
  keep <- unique(round(seq(
    1,
    dim(post_warmup)[1],
    length.out = draws_per_chain
  )))
  post_warmup <- post_warmup[keep, , , drop = FALSE]
}
retained_draws <- dim(post_warmup)[1] * dim(post_warmup)[2]

stage1 <- readRDS(stage1_file)
profile <- stage1$profiles$log_q_ats
if (is.null(profile) || !inherits(profile, "RTMBprof")) {
  stop("Stage 1 output does not contain the log_q_ats RTMBprof profile.", call. = FALSE)
}

message(
  "Evaluating ", retained_draws, " post-warmup draws from ",
  dim(post_warmup)[2], " chains."
)
mcmc_components <- RTMBprof::evaluate_mcmc_components(
  obj = obj,
  draws = post_warmup,
  components = pollock_profile_components,
  variables = parameter_names,
  objective = "joint"
)
mcmc_summary <- RTMBprof::mcmc_profile(
  mcmc_components,
  parameter = "log_q_ats",
  bins = bins,
  probs = c(0.025, 0.5, 0.975)
)

monitor <- as.data.frame(fit$monitor)
qats_monitor <- monitor[monitor$variable == "log_q_ats", , drop = FALSE]
posterior_correlations <- data.frame(
  parameter = c("log_q_avo", "log_q_cpue", "log_Rzero"),
  correlation_with_log_q_ats = vapply(
    c("log_q_avo", "log_q_cpue", "log_Rzero"),
    function(parameter) {
      stats::cor(mcmc_components$log_q_ats, mcmc_components[[parameter]])
    },
    numeric(1)
  )
)
diagnostics <- data.frame(
  chains = dim(post_warmup)[2],
  warmup_removed_per_chain = warmup,
  retained_draws = retained_draws,
  full_post_warmup_draws = total_draws,
  fitted_dimension = length(obj$par),
  gaussian_expected_objective_gap = length(obj$par) / 2,
  median_objective_gap_from_profile = stats::median(
    mcmc_components$value - min(profile$value)
  ),
  max_abs_closure_error = max(abs(mcmc_components$closure_error)),
  log_q_ats_rhat = if ("rhat" %in% names(qats_monitor)) qats_monitor$rhat else NA_real_,
  log_q_ats_ess_bulk = if ("ess_bulk" %in% names(qats_monitor)) {
    qats_monitor$ess_bulk
  } else {
    NA_real_
  },
  log_q_ats_ess_tail = if ("ess_tail" %in% names(qats_monitor)) {
    qats_monitor$ess_tail
  } else {
    NA_real_
  }
)

saveRDS(
  list(
    created = Sys.time(),
    source_mcmc_file = normalizePath(mcmc_file),
    stage1_file = normalizePath(stage1_file),
    diagnostics = diagnostics,
    posterior_correlations = posterior_correlations,
    mcmc_components = mcmc_components,
    mcmc_summary = mcmc_summary,
    profile = profile
  ),
  file.path(output_dir, "rtmbprof_stage2.rds")
)
utils::write.csv(
  diagnostics,
  file.path(output_dir, "diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  posterior_correlations,
  file.path(output_dir, "posterior_correlations.csv"),
  row.names = FALSE
)
utils::write.csv(
  mcmc_summary,
  file.path(output_dir, "log_q_ats_mcmc_summary.csv"),
  row.names = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  objective_plot <- plot(
    mcmc_components,
    profile = profile,
    summary = mcmc_summary
  ) +
    ggplot2::labs(
      x = "log acoustic-trawl survey catchability",
      y = "Change in joint negative log posterior"
    )
  ggplot2::ggsave(
    file.path(output_dir, "log_q_ats_profile_and_posterior.png"),
    objective_plot,
    width = 8,
    height = 5,
    dpi = 180
  )

  shape_plot <- RTMBprof::plot_mcmc_shape(
    mcmc_components,
    profile = profile,
    summary = mcmc_summary
  ) +
    ggplot2::labs(x = "log acoustic-trawl survey catchability")
  ggplot2::ggsave(
    file.path(output_dir, "log_q_ats_centred_shape.png"),
    shape_plot,
    width = 8,
    height = 5,
    dpi = 180
  )
}

message(
  "Wrote RTMBprof Stage 2 output to: ",
  file.path(output_dir, "rtmbprof_stage2.rds")
)
