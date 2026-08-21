# Common tidy output interface for EBS pollock model implementations.
#
# This file intentionally depends only on small, stable R interfaces. Source it
# before calling generics::tidy(), generics::glance(), or generics::augment().

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

empty_pollock_parameters <- function() {
  data.frame(
    term = character(), estimate = numeric(), std.error = numeric(),
    gradient = numeric(), estimation_type = character(), module = character(),
    year = integer(), age = integer(), engine = character(), model = character(),
    stringsAsFactors = FALSE
  )
}

empty_pollock_observations <- function() {
  data.frame(
    data_type = character(), fleet = character(), year = integer(), age = integer(),
    .truth = numeric(), .pred = numeric(), uncertainty = numeric(),
    uncertainty_scale = character(), sample_size = numeric(),
    distribution = character(), engine = character(), model = character(),
    stringsAsFactors = FALSE
  )
}

as_plain_data_frame <- function(x, template) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in setdiff(names(template), names(x))) x[[nm]] <- template[[nm]]
  x[names(template)]
}

new_pollock_fit <- function(parameters = empty_pollock_parameters(),
                            observations = empty_pollock_observations(),
                            model_summary,
                            metadata = list()) {
  stopifnot(is.list(model_summary), length(model_summary) > 0L)
  parameters <- as_plain_data_frame(parameters, empty_pollock_parameters())
  observations <- as_plain_data_frame(observations, empty_pollock_observations())
  model_summary <- as.data.frame(model_summary, stringsAsFactors = FALSE)
  if (nrow(model_summary) != 1L) stop("model_summary must contain exactly one row.")
  structure(
    list(
      parameters = parameters,
      observations = observations,
      model_summary = model_summary,
      metadata = metadata
    ),
    class = "pollock_fit"
  )
}

tidy.pollock_fit <- function(x,
                             parameters = c("fixed_effect", "random_effect"),
                             conf.int = FALSE,
                             conf.level = 0.95,
                             ...) {
  valid <- c("fixed_effect", "random_effect", "derived_quantity", "fixed_input")
  bad <- setdiff(parameters, valid)
  if (length(bad)) stop("Unknown estimation type: ", paste(bad, collapse = ", "))
  out <- x$parameters[x$parameters$estimation_type %in% parameters, , drop = FALSE]
  out$statistic <- out$estimate / out$std.error
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  if (conf.int) {
    z <- stats::qnorm((1 + conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }
  tibble::as_tibble(out)
}

glance.pollock_fit <- function(x, ...) tibble::as_tibble(x$model_summary)

augment.pollock_fit <- function(x, include_weights = FALSE, ...) {
  out <- x$observations
  out <- out[is.finite(out$.truth) & is.finite(out$.pred), , drop = FALSE]
  if (include_weights) {
    out$.weight <- ifelse(
      out$uncertainty_scale == "observation" &
        is.finite(out$uncertainty) & out$uncertainty > 0,
      1 / out$uncertainty^2,
      NA_real_
    )
  }
  tibble::as_tibble(out)
}

pollock_fit_metrics <- function(
    x,
    group_by = c("engine", "model", "data_type", "fleet"),
    data_types = c("catch", "index"),
    metrics = NULL,
    weighted = FALSE) {
  if (!requireNamespace("yardstick", quietly = TRUE)) {
    stop("Package 'yardstick' is required for pollock_fit_metrics().")
  }
  aug <- augment.pollock_fit(x, include_weights = weighted)
  aug <- aug[
    aug$data_type %in% data_types & is.finite(aug$.truth) & is.finite(aug$.pred),
    , drop = FALSE
  ]
  if (!nrow(aug)) return(tibble::tibble())
  bad <- setdiff(group_by, names(aug))
  if (length(bad)) stop("Unknown grouping column: ", paste(bad, collapse = ", "))
  if (is.null(metrics)) {
    metrics <- yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq_trad)
  }
  aug <- dplyr::group_by(aug, dplyr::across(dplyr::all_of(group_by)))
  if (weighted) {
    keep <- is.finite(aug$.weight) & aug$.weight > 0
    aug <- aug[keep, , drop = FALSE]
    if (!nrow(aug)) stop("No finite inverse-variance weights are available.")
    aug$.weight <- hardhat::importance_weights(aug$.weight)
    metrics(aug, truth = .truth, estimate = .pred, case_weights = .weight)
  } else {
    metrics(aug, truth = .truth, estimate = .pred)
  }
}

indexed_parameter_rows <- function(x, engine, model, estimation_type = "fixed_effect") {
  if (is.null(x) || !length(x)) return(empty_pollock_parameters())
  rows <- lapply(names(x), function(term) {
    value <- x[[term]]
    if (!is.numeric(value) || !length(value)) return(NULL)
    idx <- if (length(value) == 1L) "" else paste0("[", seq_along(value), "]")
    data.frame(
      term = paste0(term, idx), estimate = as.numeric(value), std.error = NA_real_,
      gradient = NA_real_, estimation_type = estimation_type, module = NA_character_,
      year = NA_integer_, age = NA_integer_, engine = engine, model = model,
      stringsAsFactors = FALSE
    )
  })
  as_plain_data_frame(do.call(rbind, rows), empty_pollock_parameters())
}

derived_time_series_rows <- function(values, term, years, engine, model) {
  if (is.null(values) || !length(values)) return(empty_pollock_parameters())
  data.frame(
    term = term, estimate = as.numeric(values), std.error = NA_real_, gradient = NA_real_,
    estimation_type = "derived_quantity", module = "population", year = as.integer(years),
    age = NA_integer_, engine = engine, model = model, stringsAsFactors = FALSE
  ) |> as_plain_data_frame(empty_pollock_parameters())
}

observation_rows <- function(truth, pred, data_type, fleet, years,
                             ages = NA_integer_, uncertainty = NA_real_,
                             uncertainty_scale = NA_character_,
                             sample_size = NA_real_, distribution = NA_character_,
                             engine, model) {
  n <- length(truth)
  if (length(pred) != n) stop("Observed and predicted vectors have different lengths.")
  recycle <- function(x) if (length(x) == 1L) rep(x, n) else x
  data.frame(
    data_type = data_type, fleet = fleet, year = recycle(years), age = recycle(ages),
    .truth = as.numeric(truth), .pred = as.numeric(pred),
    uncertainty = recycle(uncertainty),
    uncertainty_scale = recycle(uncertainty_scale), sample_size = recycle(sample_size),
    distribution = distribution, engine = engine, model = model,
    stringsAsFactors = FALSE
  ) |> as_plain_data_frame(empty_pollock_observations())
}

composition_rows <- function(obs, pred, years, sample_size, fleet, engine, model) {
  if (is.null(obs) || is.null(pred) || !length(obs)) return(empty_pollock_observations())
  obs <- as.matrix(obs); pred <- as.matrix(pred)
  if (!identical(dim(obs), dim(pred))) stop("Composition dimensions do not match.")
  observation_rows(
    truth = as.vector(t(obs)), pred = as.vector(t(pred)),
    data_type = "age_composition", fleet = fleet,
    years = rep(as.integer(years), each = ncol(obs)),
    ages = rep(seq_len(ncol(obs)), times = nrow(obs)),
    sample_size = rep(as.numeric(sample_size), each = ncol(obs)),
    distribution = "composition", engine = engine, model = model
  )
}

as_pollock_fit_rtmb <- function(x, model = x$label %||% "base") {
  report <- x$report %||% x$rtmb
  if (is.null(report)) stop("RTMB saved object has no report component.")
  engine <- "custom RTMB"
  years <- suppressWarnings(as.integer(rownames(report$N)))
  if (!length(years) || anyNA(years)) years <- seq_along(report$obs_catch)
  if (!is.null(x$parameter_table)) {
    fixed <- as_plain_data_frame(x$parameter_table, empty_pollock_parameters())
    fixed$engine <- engine
    fixed$model <- model
  } else {
    # Legacy files retain a full parameter state but not the active RTMB map.
    # Treat these values as inputs rather than falsely labelling every element
    # as an estimated fixed effect.
    fixed <- indexed_parameter_rows(
      x$fixed_parameters %||% list(), engine, model,
      estimation_type = "fixed_input"
    )
  }
  derived <- rbind(
    derived_time_series_rows(report$SSB, "spawning_biomass", years, engine, model),
    derived_time_series_rows(report$N[, 1], "recruitment", years, engine, model)
  )
  obs <- rbind(
    observation_rows(report$obs_catch, report$pred_catch, "catch", "fishery", years,
                     distribution = "lognormal", engine = engine, model = model),
    composition_rows(report$oac_fsh, report$phat_fsh, report$yrs_fsh_data,
                     report$sam_fsh, "fishery", engine, model),
    composition_rows(report$oac_bts, report$phat_bts, report$yrs_bts_data,
                     report$sam_bts, "BTS", engine, model),
    composition_rows(report$oac_ats, report$phat_ats, report$yrs_ats_data,
                     report$sam_ats, "ATS", engine, model)
  )
  n_random <- as.integer(x$n_random %||% 0L)
  objective <- as.numeric(x$objective %||% report$tot_like %||% NA_real_)
  joint_nll <- as.numeric(x$evaluated_total %||% report$tot_like %||% NA_real_)
  summary <- list(
    engine = engine, model = model,
    objective_type = if (n_random > 0L) "marginal" else "joint",
    marginal_nll = if (n_random > 0L) objective else NA_real_,
    joint_nll = joint_nll, npar_fixed = as.integer((x$n_parameters %||% nrow(fixed)) - n_random),
    npar_random = n_random, max_gradient = as.numeric(x$max_gradient %||% NA_real_),
    convergence = as.integer(x$convergence %||% NA_integer_),
    converged = isTRUE((x$convergence %||% 1L) == 0L) &&
      isTRUE((x$max_gradient %||% Inf) < 0.001),
    hessian_positive_definite = x$hessian_positive_definite %||% NA,
    terminal_ssb = tail(as.numeric(report$SSB), 1), runtime_secs = as.numeric(x$seconds %||% NA_real_)
  )
  new_pollock_fit(rbind(fixed, derived), obs, summary, metadata = x[c("form", "label")])
}

as_pollock_fit_rceattle <- function(x, model = "nonparametric_pm") {
  if (inherits(x, "Rceattle")) fit <- x else fit <- x[[model]] %||% x
  if (is.null(fit$quantities) || is.null(fit$data_list)) stop("Unrecognized Rceattle fit.")
  engine <- "Rceattle"
  diagnostics <- fit$opt$diagnostics
  se <- if (!is.null(fit$sdrep$cov.fixed)) sqrt(pmax(base::diag(fit$sdrep$cov.fixed), 0)) else NA_real_
  fixed <- data.frame(
    term = make.unique(as.character(diagnostics$Param)), estimate = as.numeric(diagnostics$MLE),
    std.error = se, gradient = as.numeric(diagnostics$final_gradient),
    estimation_type = "fixed_effect", module = NA_character_, year = NA_integer_,
    age = NA_integer_, engine = engine, model = model, stringsAsFactors = FALSE
  ) |> as_plain_data_frame(empty_pollock_parameters())
  years <- fit$data_list$styr:fit$data_list$projyr
  derived <- rbind(
    derived_time_series_rows(fit$quantities$ssb[1, ], "spawning_biomass", years, engine, model),
    derived_time_series_rows(fit$quantities$R[1, ], "recruitment", years, engine, model)
  )
  catch <- fit$data_list$catch_data
  index <- fit$data_list$index_data
  obs <- rbind(
    observation_rows(catch$Catch, fit$quantities$catch_hat, "catch", catch$Fleet_name,
                     catch$Year, uncertainty = catch$Log_sd, uncertainty_scale = "log",
                     distribution = "lognormal", engine = engine, model = model),
    observation_rows(index$Observation, fit$quantities$index_hat, "index", index$Fleet_name,
                     index$Year, uncertainty = index$Log_sd, uncertainty_scale = "log",
                     distribution = "lognormal", engine = engine, model = model)
  )
  opt <- fit$opt
  summary <- list(
    engine = engine, model = model, objective_type = "joint",
    marginal_nll = NA_real_, joint_nll = as.numeric(fit$quantities$jnll),
    npar_fixed = as.integer(opt$number_of_coefficients[["Fixed"]]),
    npar_random = as.integer(opt$number_of_coefficients[["Random"]]),
    max_gradient = as.numeric(opt$max_gradient), convergence = NA_integer_,
    converged = isTRUE(opt$max_gradient < 0.001),
    hessian_positive_definite = fit$sdrep$pdHess %||% NA,
    terminal_ssb = tail(as.numeric(fit$quantities$ssb), 1),
    runtime_secs = as.numeric(opt$time_for_run)
  )
  new_pollock_fit(rbind(fixed, derived), obs, summary)
}

as_pollock_fit_sporc <- function(x, model = x$case$id %||% "base") {
  if (is.null(x$rep) || is.null(x$data)) stop("Unrecognized SPoRC fit.")
  engine <- "SPoRC"
  par <- x$sdrep$par.fixed %||% x$optim$par
  se <- if (!is.null(x$sdrep$cov.fixed)) sqrt(pmax(diag(x$sdrep$cov.fixed), 0)) else NA_real_
  grad <- as.numeric(x$sdrep$gradient.fixed %||% NA_real_)
  fixed <- data.frame(
    term = make.unique(names(par)), estimate = as.numeric(par), std.error = se,
    gradient = rep(grad, length.out = length(par)), estimation_type = "fixed_effect",
    module = NA_character_, year = NA_integer_, age = NA_integer_, engine = engine,
    model = model, stringsAsFactors = FALSE
  ) |> as_plain_data_frame(empty_pollock_parameters())
  years <- as.integer(x$data$years)
  derived <- rbind(
    derived_time_series_rows(as.numeric(x$rep$SSB), "spawning_biomass", years, engine, model),
    derived_time_series_rows(as.numeric(x$rep$Rec), "recruitment", years, engine, model)
  )
  make_array_stream <- function(obs, pred, data_type, prefix, uncertainty = NA_real_) {
    obs <- as.array(obs); pred <- as.array(pred)
    to_year_matrix <- function(value) {
      d <- dim(value)
      year_dim <- which(d == length(years))[1]
      if (is.na(year_dim)) stop(prefix, " has no dimension matching model years.")
      perm <- c(year_dim, setdiff(seq_along(d), year_dim))
      matrix(aperm(value, perm), nrow = length(years))
    }
    o <- to_year_matrix(obs); p <- to_year_matrix(pred)
    if (!identical(dim(o), dim(p))) stop(prefix, " dimensions do not match after removing singleton dimensions.")
    n_fleet <- ncol(o)
    uncertainty <- if (length(uncertainty) == 1L) {
      uncertainty
    } else {
      as.vector(to_year_matrix(as.array(uncertainty)))
    }
    observation_rows(
      as.vector(o), as.vector(p), data_type,
      rep(paste(prefix, seq_len(n_fleet)), each = length(years)),
      rep(years, times = n_fleet), uncertainty = as.vector(uncertainty),
      uncertainty_scale = if (length(uncertainty) == 1L && is.na(uncertainty)) NA_character_ else "log",
      distribution = "lognormal", engine = engine, model = model
    )
  }
  obs <- rbind(
    make_array_stream(x$data$ObsCatch, x$rep$PredCatch, "catch", "fishery"),
    make_array_stream(x$data$ObsFishIdx, x$rep$PredFishIdx, "index", "fishery index",
                      x$data$ObsFishIdx_SE),
    make_array_stream(x$data$ObsSrvIdx, x$rep$PredSrvIdx, "index", "survey",
                      x$data$ObsSrvIdx_SE)
  )
  opt <- x$optim
  n_random <- length(x$env$random %||% integer())
  summary <- list(
    engine = engine, model = model,
    objective_type = if (n_random > 0L) "marginal" else "joint",
    marginal_nll = if (n_random > 0L) as.numeric(opt$objective) else NA_real_,
    joint_nll = if (n_random > 0L) NA_real_ else as.numeric(opt$objective),
    npar_fixed = length(par), npar_random = n_random,
    max_gradient = max(abs(grad), na.rm = TRUE), convergence = as.integer(opt$convergence),
    converged = isTRUE(opt$convergence == 0L) && isTRUE(max(abs(grad), na.rm = TRUE) < 0.001),
    hessian_positive_definite = x$sdrep$pdHess %||% NA,
    terminal_ssb = tail(as.numeric(x$rep$SSB), 1), runtime_secs = NA_real_
  )
  new_pollock_fit(rbind(fixed, derived), obs, summary, metadata = x$case %||% list())
}
