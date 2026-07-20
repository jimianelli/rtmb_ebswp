# Likelihood profiles by reported objective component.
#
# These helpers adapt the SBT component-profile workflow to the EBS pollock
# RTMB model. They assume the objective reports likelihood pieces via REPORT()
# in R/Rpm.R, and work with an RTMB MakeADFun object.

profile_get_report_value <- function(report, names, default = NA_real_) {
  for (nm in names) {
    if (!is.null(report[[nm]])) {
      return(report[[nm]])
    }
  }
  default
}

profile_component_sum <- function(report, names) {
  value <- profile_get_report_value(report, names)
  if (is.null(value) || all(is.na(value))) {
    return(NA_real_)
  }
  sum(value, na.rm = TRUE)
}

pollock_profile_components <- function(report) {
  age_like <- profile_get_report_value(report, "age_like")
  age_total <- attr(age_like, "total")
  if (is.null(age_total)) {
    age_total <- sum(age_like, na.rm = TRUE)
  }

  c(
    Prior = profile_component_sum(report, "Priors"),
    Recruitment = profile_component_sum(
      report,
      c("rec_like_component", "rec_like$rec_like", "rec_like.rec_like", "rec_like")
    ),
    Age = age_total,
    BTS = profile_component_sum(report, "bts_like"),
    ATS = profile_component_sum(report, "ats_like"),
    `ATS age1` = profile_component_sum(report, "ats_age1_like"),
    CPUE = profile_component_sum(report, "cpue_like"),
    AVO = profile_component_sum(report, "avo_like"),
    Catch = profile_component_sum(report, "cat_like"),
    Fpen = profile_component_sum(report, "Fpen_like"),
    Selectivity = profile_component_sum(report, "sel_like"),
    `Selectivity dev` = profile_component_sum(report, "sel_like_dev"),
    Weight = profile_component_sum(report, "wt_like")
  )
}

profile_start_par <- function(obj, par = NULL) {
  if (!is.null(par)) {
    return(par)
  }
  if (!is.null(obj$env$last.par.best)) {
    par <- obj$env$last.par.best
  } else {
    par <- obj$par
  }
  if (!is.null(obj$env$random)) {
    par <- par[-obj$env$random]
  }
  par
}

profile_resolve_parameter <- function(par, name) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be one non-empty parameter name.", call. = FALSE)
  }

  indexed <- regexec("^(.+)\\[([0-9]+)\\]$", name)
  parts <- regmatches(name, indexed)[[1]]
  if (length(parts) > 0L) {
    base_name <- parts[2]
    occurrence <- as.integer(parts[3])
    matches <- which(names(par) == base_name)
    if (occurrence < 1L || occurrence > length(matches)) {
      stop(
        "Parameter occurrence not found: ", name,
        ". Available occurrences: ", length(matches),
        call. = FALSE
      )
    }
    return(list(index = matches[occurrence], label = name, base_name = base_name))
  }

  matches <- which(names(par) == name)
  if (length(matches) == 0L) {
    stop("Parameter not found: ", name, call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop(
      "Parameter name is repeated: ", name, ". Use ", name,
      "[n], where n is 1 through ", length(matches), ".",
      call. = FALSE
    )
  }
  list(index = matches, label = name, base_name = name)
}

profile_validate_values <- function(values) {
  if (!is.numeric(values) || length(values) < 1L || any(!is.finite(values))) {
    stop("`values` must contain one or more finite numeric values.", call. = FALSE)
  }
  as.numeric(values)
}

profile_result_row <- function(parameter, parameter_index, value, objective,
                               report, component_fun, ...) {
  components <- component_fun(report)
  component_total <- sum(components, na.rm = TRUE)
  data.frame(
    parameter = parameter,
    parameter_index = parameter_index,
    value = value,
    objective = objective,
    ...,
    as.list(components),
    Other = objective - component_total,
    check.names = FALSE
  )
}

profile_components_slice <- function(obj, name, values,
                                     component_fun = pollock_profile_components,
                                     par = NULL) {
  par0 <- profile_start_par(obj, par = par)
  target <- profile_resolve_parameter(par0, name)
  idx <- target$index
  values <- profile_validate_values(values)
  on.exit(obj$fn(par0), add = TRUE)

  out <- vector("list", length(values))
  for (i in seq_along(values)) {
    par_i <- par0
    par_i[idx] <- values[i]
    objective_i <- obj$fn(par_i)
    report_i <- obj$report(par_i)

    out[[i]] <- profile_result_row(
      parameter = target$label,
      parameter_index = idx,
      value = values[i],
      objective = objective_i,
      report = report_i,
      component_fun = component_fun
    )
  }

  do.call(rbind, out)
}

profile_components_reopt <- function(obj, name, values,
                                     component_fun = pollock_profile_components,
                                     par = NULL,
                                     control = list(iter.max = 500, eval.max = 500),
                                     gradient_tolerance = NULL,
                                     max_polish = 3L,
                                     trace = TRUE) {
  par0 <- profile_start_par(obj, par = par)
  target <- profile_resolve_parameter(par0, name)
  idx <- target$index
  values <- profile_validate_values(values)
  on.exit(obj$fn(par0), add = TRUE)

  free_idx <- setdiff(seq_along(par0), idx)
  if (length(free_idx) == 0L) {
    stop("No free parameters remain after fixing `", name, "`.", call. = FALSE)
  }
  start_free <- par0[free_idx]
  out <- vector("list", length(values))

  for (i in seq_along(values)) {
    fixed_value <- values[i]

    fn_free <- function(free_par) {
      par_i <- par0
      par_i[idx] <- fixed_value
      par_i[free_idx] <- free_par
      obj$fn(par_i)
    }

    gr_free <- function(free_par) {
      par_i <- par0
      par_i[idx] <- fixed_value
      par_i[free_idx] <- free_par
      obj$gr(par_i)[free_idx]
    }

    opt_i <- nlminb(
      start = start_free,
      objective = fn_free,
      gradient = gr_free,
      control = control
    )

    max_gradient <- max(abs(gr_free(opt_i$par)))
    polish_steps <- 0L
    while (!is.null(gradient_tolerance) &&
           is.finite(max_gradient) &&
           max_gradient > gradient_tolerance &&
           polish_steps < max_polish) {
      previous_gradient <- max_gradient
      polished <- nlminb(
        start = opt_i$par,
        objective = fn_free,
        gradient = gr_free,
        control = control
      )
      polished_gradient <- max(abs(gr_free(polished$par)))
      if (polished$objective <= opt_i$objective + 1e-8) {
        opt_i <- polished
        max_gradient <- polished_gradient
      }
      polish_steps <- polish_steps + 1L
      if (max_gradient >= previous_gradient) {
        break
      }
    }

    par_i <- par0
    par_i[idx] <- fixed_value
    par_i[free_idx] <- opt_i$par
    report_i <- obj$report(par_i)

    if (isTRUE(trace)) {
      message(
        "Profile ", name, " = ", signif(fixed_value, 6),
        "; objective = ", signif(opt_i$objective, 8),
        "; convergence = ", opt_i$convergence
      )
    }

    out[[i]] <- profile_result_row(
      parameter = target$label,
      parameter_index = idx,
      value = fixed_value,
      objective = opt_i$objective,
      report = report_i,
      component_fun = component_fun,
      convergence = opt_i$convergence,
      evaluations = unname(opt_i$evaluations[["function"]]),
      max_gradient = max_gradient,
      polish_steps = polish_steps
    )

    start_free <- opt_i$par
  }

  do.call(rbind, out)
}

profile_components_long <- function(profile, rescale = TRUE,
                                    include_objective = TRUE) {
  reserved <- intersect(
    c(
      "parameter", "parameter_index", "value", "objective", "convergence",
      "evaluations", "max_gradient", "polish_steps"
    ),
    names(profile)
  )
  component_data <- profile[setdiff(names(profile), reserved)]
  if (isTRUE(include_objective)) {
    component_data <- cbind(
      data.frame(Objective = profile$objective),
      component_data
    )
  }
  component_cols <- names(component_data)

  out <- data.frame(
    profile[rep(seq_len(nrow(profile)), length(component_cols)), reserved, drop = FALSE],
    component = rep(component_cols, each = nrow(profile)),
    nll = as.numeric(unlist(component_data, use.names = FALSE)),
    row.names = NULL,
    check.names = FALSE
  )
  out$component <- factor(out$component, levels = component_cols)

  if (isTRUE(rescale)) {
    split_index <- split(seq_len(nrow(out)), out$component)
    out$delta <- out$nll
    for (idx in split_index) {
      out$delta[idx] <- out$nll[idx] - min(out$nll[idx], na.rm = TRUE)
    }
  }

  out
}

plot_profile_components <- function(profile, rescale = TRUE,
                                    y_limits = c(0, 2.1),
                                    likelihood_threshold = 1.92) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }
  if (!requireNamespace("ggthemes", quietly = TRUE)) {
    stop("Package `ggthemes` is required for plotting.", call. = FALSE)
  }
  if (!is.numeric(y_limits) || length(y_limits) != 2L ||
      any(!is.finite(y_limits)) || y_limits[1] >= y_limits[2]) {
    stop("`y_limits` must contain two increasing finite values.", call. = FALSE)
  }

  dat <- profile_components_long(
    profile, rescale = rescale, include_objective = TRUE
  )
  yvar <- if (isTRUE(rescale)) "delta" else "nll"
  ylab <- if (isTRUE(rescale)) "Delta NLL" else "NLL"
  xlab <- unique(profile$parameter)
  if (length(xlab) != 1L) {
    xlab <- "Profile value"
  }
  fitted_value <- profile$value[which.min(profile$objective)]
  threshold_data <- data.frame(
    component = factor("Objective", levels = levels(dat$component)),
    threshold = likelihood_threshold
  )

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$value, y = .data[[yvar]])) +
    ggplot2::geom_hline(
      data = threshold_data,
      ggplot2::aes(yintercept = .data$threshold),
      linetype = 2,
      color = "grey35",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = fitted_value,
      linetype = 3,
      color = "grey45"
    ) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(ggplot2::vars(.data$component), scales = "fixed") +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      caption = paste0(
        "Vertical dotted line: lowest evaluated objective; ",
        "horizontal dashed line: Delta NLL = ", likelihood_threshold,
        " (Objective panel only)."
      )
    ) +
    ggthemes::theme_few()
}
