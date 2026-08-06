# Run retrospective fits for the RTMB-ADMB bridge model.
# Usage: Rscript R/run_retrospective.R

source(file.path("R", "config.R"))

rtmb_dir <- normalizePath(getwd(), mustWork = TRUE)
base_file <- file.path(rtmb_dir, "analysis", "output", "base.rds")
if (!file.exists(base_file)) {
  stop("Run R/write_output.R before the retrospective analysis.")
}
base_saved <- readRDS(base_file)
if (!identical(
  base_saved$metadata$configuration,
  "September 2025 fixed-parameter ADMB-to-RTMB bridge"
)) {
  stop("base.rds is not identified as the September 2025 bridge configuration.")
}
base_md5 <- unname(tools::md5sum(base_file))
fixed_params <- names(map_obj)[vapply(
  map_obj, function(x) all(is.na(x)), logical(1)
)]

trim_rows <- function(x, keep) {
  if (is.null(x)) {
    return(x)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    return(x[keep, , drop = FALSE])
  }
  if (is.array(x) && length(dim(x)) >= 2) {
    idx <- rep(list(TRUE), length(dim(x)))
    idx[[1]] <- keep
    return(do.call(`[`, c(list(x), idx, list(drop = FALSE))))
  }
  x[keep]
}

trim_rtmb_data <- function(data, peel) {
  out <- data
  full_endyr <- data$endyr
  new_endyr <- full_endyr - peel
  if (new_endyr <= data$styr) {
    stop("Peel leaves too few years: ", peel)
  }

  year_keep <- data$styr:new_endyr
  year_index <- seq_along(data$styr:full_endyr)[year_keep - data$styr + 1]

  out[names(out) == "endyr"] <- rep(list(new_endyr), sum(names(out) == "endyr"))
  out[names(out) == "endyr_est"] <- rep(list(new_endyr - out$omitSR), sum(names(out) == "endyr_est"))
  out$nyrs <- length(year_keep)
  out$endyr_wt <- min(out$endyr_wt, new_endyr)

  annual_names <- c(
    "wt_fsh", "wt_ssb", "obs_catch", "obs_effort", "ew_ind"
  )
  for (nm in annual_names) {
    if (is.null(out[[nm]])) {
      next
    }
    n_available <- if (is.null(dim(out[[nm]]))) length(out[[nm]]) else dim(out[[nm]])[1]
    if (n_available >= max(year_index)) {
      out[[nm]] <- trim_rows(out[[nm]], year_index)
    }
  }

  keep_cpue <- out$yrs_cpue <= new_endyr
  out$yrs_cpue <- out$yrs_cpue[keep_cpue]
  out$obs_cpue <- out$obs_cpue[keep_cpue]
  out$obs_cpue_std <- out$obs_cpue_std[keep_cpue]
  out$obs_cpue_var <- out$obs_cpue_std^2
  out$n_cpue <- length(out$yrs_cpue)

  keep_avo <- out$yrs_avo <= new_endyr
  out$yrs_avo <- out$yrs_avo[keep_avo]
  out$ob_avo <- out$ob_avo[keep_avo]
  out$ob_avo_std <- out$ob_avo_std[keep_avo]
  out$obs_avo_var <- out$ob_avo_std^2
  out$wt_avo <- trim_rows(out$wt_avo, keep_avo)
  out$n_avo <- length(out$yrs_avo)

  keep_fsh <- out$yrs_fsh_data <= new_endyr
  out$yrs_fsh_data <- out$yrs_fsh_data[keep_fsh]
  out$sam_fsh <- out$sam_fsh[keep_fsh]
  out$err_fsh <- out$err_fsh[keep_fsh]
  out$oac_fsh_data <- trim_rows(out$oac_fsh_data, keep_fsh)
  out$oac_fsh <- trim_rows(out$oac_fsh, keep_fsh)
  out$n_fsh <- length(out$yrs_fsh_data)

  keep_bts <- out$yrs_bts_data <= new_endyr
  out$yrs_bts_data <- out$yrs_bts_data[keep_bts]
  out$sam_bts <- out$sam_bts[keep_bts]
  out$err_bts <- out$err_bts[keep_bts]
  out$ob_bts <- out$ob_bts[keep_bts]
  out$ob_bts_std <- out$ob_bts_std[keep_bts]
  out$std_ob_bts <- out$std_ob_bts[keep_bts]
  out$ot_bts <- out$ot_bts[keep_bts]
  out$wt_bts <- trim_rows(out$wt_bts, keep_bts)
  out$oac_bts <- trim_rows(out$oac_bts, keep_bts)
  out$cov_matrix <- out$cov_matrix[keep_bts, keep_bts, drop = FALSE]
  out$inv_bts_cov <- solve(out$cov_matrix)
  out$n_bts <- length(out$yrs_bts_data)

  keep_ats <- out$yrs_ats_data <= new_endyr
  out$yrs_ats_data <- out$yrs_ats_data[keep_ats]
  out$sam_ats <- out$sam_ats[keep_ats]
  out$err_ats <- out$err_ats[keep_ats]
  out$ob_ats <- out$ob_ats[keep_ats]
  out$ob_ats_std <- out$ob_ats_std[keep_ats]
  out$ot_ats <- out$ot_ats[keep_ats]
  out$ot_ats_std <- out$ot_ats_std[keep_ats]
  out$oa1_ats <- out$oa1_ats[keep_ats]
  out$lvarb_ats <- out$lvarb_ats[keep_ats]
  out$wt_ats <- trim_rows(out$wt_ats, keep_ats)
  out$oac_ats <- trim_rows(out$oac_ats, keep_ats)
  out$n_ats <- length(out$yrs_ats_data)

  out$nagecomp <- c(out$n_fsh, out$n_bts, out$n_ats)
  out$wt_fut <- out$wt_fsh[nrow(out$wt_fsh), ]

  if (!is.null(out$yrs_data) && !is.null(out$nyrs_data)) {
    for (i in seq_len(nrow(out$yrs_data))) {
      keep <- !is.na(out$yrs_data[i, ]) & out$yrs_data[i, ] <= new_endyr
      out$nyrs_data[i] <- sum(keep)
      out$yrs_data[i, ] <- NA
      if (out$nyrs_data[i] > 0) {
        out$yrs_data[i, seq_len(out$nyrs_data[i])] <- data$yrs_data[i, keep]
      }
    }
  }

  out
}

trim_rows_n <- function(x, n) {
  if (is.null(x)) {
    return(x)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    return(x[seq_len(n), , drop = FALSE])
  }
  if (is.array(x) && length(dim(x)) >= 2) {
    idx <- rep(list(TRUE), length(dim(x)))
    idx[[1]] <- seq_len(n)
    return(do.call(`[`, c(list(x), idx, list(drop = FALSE))))
  }
  x[seq_len(n)]
}

trim_rtmb_parms <- function(parms, data) {
  out <- parms
  nyrs <- data$endyr - data$styr + 1
  nyrs_bts <- data$endyr - data$styr_bts + 1
  nyrs_ats_dev <- length(1995:min(2024, data$endyr))
  nyrs_fsh_dev <- length(1965:min(2023, data$endyr))

  for (nm in c("log_rec_devs", "log_F_devs", "sel_dif1_fsh_dev",
               "sel_a501_fsh_dev", "sel_trm2_fsh_dev")) {
    out[[nm]] <- trim_rows_n(out[[nm]], nyrs)
  }

  out$M_dev <- trim_rows_n(out$M_dev, nyrs - 1)
  out$yr_eff <- trim_rows_n(
    out$yr_eff, data$endyr_wt - data$styr_wt + 1
  )

  for (nm in c("sel_slp_bts_dev", "sel_a50_bts_dev", "sel_age_one_bts_dev")) {
    out[[nm]] <- trim_rows_n(out[[nm]], nyrs_bts)
  }

  out$sel_devs_bts <- trim_rows_n(out$sel_devs_bts, max(nyrs_bts - 1, 1))
  out$sel_devs_ats <- trim_rows_n(out$sel_devs_ats, nyrs_ats_dev)
  out$sel_devs_fsh <- trim_rows_n(out$sel_devs_fsh, nyrs_fsh_dev)

  out
}

fit_one_peel <- function(peel) {
  message("Running RTMB retrospective peel ", peel)
  data <<- trim_rtmb_data(full_data, peel)
  parms_peel <- trim_rtmb_parms(start_parms, data)
  map_peel <- create_map_from_par(
    parms_peel, parms_peel,
    exact_names = fixed_params,
    exclude_patterns = "xxx"
  )
  obj_peel <- RTMB::MakeADFun(rpm, parms_peel, map = map_peel, silent = TRUE)
  fit <- nlminb(
    obj_peel$par,
    obj_peel$fn,
    obj_peel$gr,
    control = list(eval.max = 5000, iter.max = 3000)
  )
  first_gradient <- max(abs(obj_peel$gr(fit$par)), na.rm = TRUE)
  optimization_passes <- 1L
  if (is.finite(first_gradient) && first_gradient > 1e-3) {
    fit_restart <- nlminb(
      fit$par,
      obj_peel$fn,
      obj_peel$gr,
      control = list(
        eval.max = 5000,
        iter.max = 3000,
        rel.tol = 1e-12,
        x.tol = 1e-10
      )
    )
    if (is.finite(fit_restart$objective) &&
        fit_restart$objective <= fit$objective + 1e-8) {
      fit <- fit_restart
      optimization_passes <- 2L
    }
  }
  obj_peel$fn(fit$par)
  report <- obj_peel$report()
  fitted_parms <- obj_peel$env$parList(fit$par)
  start_parms <<- fitted_parms
  if (!is.list(report) || is.null(report$SSB) || is.null(report$N)) {
    data_report <- data
    data_report$return_nll_only <- 0
    data_report[names(data_report) == "return_nll_only"] <-
      rep(list(0), sum(names(data_report) == "return_nll_only"))
    data <<- data_report
    report <- rpm(fitted_parms)
    data_report$return_nll_only <- 1
    data_report[names(data_report) == "return_nll_only"] <-
      rep(list(1), sum(names(data_report) == "return_nll_only"))
    data <<- data_report
  }
  if (is.list(report) && !is.null(report$rtmb)) {
    report <- report$rtmb
  }

  list(
    peel = peel,
    terminal_year = data$endyr,
    convergence = fit$convergence,
    objective = fit$objective,
    max_gradient = max(abs(obj_peel$gr(fit$par)), na.rm = TRUE),
    optimization_passes = optimization_passes,
    report = report
  )
}

full_data <- data
full_parms <- parms
start_parms <- full_parms
peel_spec <- Sys.getenv("RTMB_RETRO_PEELS", unset = "0:9")
peels <- eval(parse(text = peel_spec))
stop_on_error <- tolower(Sys.getenv("RTMB_RETRO_STOP_ON_ERROR", "false")) %in%
  c("1", "true", "yes")
retro <- lapply(peels, function(peel) {
  if (stop_on_error) {
    fit_one_peel(peel)
  } else {
    tryCatch(
      fit_one_peel(peel),
      error = function(e) {
      message("RTMB retrospective peel ", peel, " failed: ", conditionMessage(e))
      list(
        peel = peel,
        terminal_year = full_data$endyr - peel,
        convergence = NA_integer_,
        objective = NA_real_,
        max_gradient = NA_real_,
        error = conditionMessage(e),
        report = NULL
      )
      }
    )
  }
})

series <- do.call(rbind, lapply(retro, function(x) {
  if (is.null(x$report) || !is.list(x$report)) {
    return(NULL)
  }
  yrs <- full_data$styr:x$terminal_year
  if (length(x$report$SSB) != length(yrs) || nrow(x$report$N) != length(yrs)) {
    return(NULL)
  }
  data.frame(
    peel = x$peel,
    terminal_year = x$terminal_year,
    year = yrs,
    SSB = as.numeric(x$report$SSB),
    Recruitment = as.numeric(x$report$N[, 1]),
    stringsAsFactors = FALSE
  )
}))

diagnostics <- do.call(rbind, lapply(retro, function(x) {
  data.frame(
    peel = x$peel,
    terminal_year = x$terminal_year,
    convergence = x$convergence,
    objective = x$objective,
    max_gradient = x$max_gradient,
    optimization_passes = if (!is.null(x$optimization_passes)) {
      x$optimization_passes
    } else NA_integer_,
    error = if (!is.null(x$error)) x$error else NA_character_,
    stringsAsFactors = FALSE
  )
}))

mohn <- NULL
if (!is.null(series) && nrow(series) > 0 && any(series$peel == 0)) {
  full_series <- series[series$peel == 0, c("year", "SSB", "Recruitment")]
  peel_terminal <- series[series$peel > 0 & series$year == series$terminal_year, ]
  comparisons <- merge(
    peel_terminal,
    full_series,
    by = "year",
    suffixes = c("_peel", "_full")
  )
  comparisons$SSB_relative_difference <-
    (comparisons$SSB_peel - comparisons$SSB_full) / comparisons$SSB_full
  comparisons$Recruitment_relative_difference <-
    (comparisons$Recruitment_peel - comparisons$Recruitment_full) /
    comparisons$Recruitment_full
  mohn <- list(
    comparisons = comparisons,
    rho = data.frame(
      quantity = c("SSB", "Recruitment"),
      rho = c(
        mean(comparisons$SSB_relative_difference, na.rm = TRUE),
        mean(comparisons$Recruitment_relative_difference, na.rm = TRUE)
      ),
      n_peels = nrow(comparisons),
      stringsAsFactors = FALSE
    )
  )
}

out <- list(
  created = Sys.time(),
  base_lineage = list(
    file = normalizePath(base_file),
    md5 = base_md5,
    created = base_saved$metadata$created,
    configuration = base_saved$metadata$configuration,
    bridge_metrics = base_saved$metadata$bridge_metrics
  ),
  peels = peels,
  series = series,
  diagnostics = diagnostics,
  mohn = mohn,
  fits = retro
)

max_peel <- max(peels)
output_path <- file.path(
  rtmb_dir, "analysis", "output", sprintf("retro_%d_peel.rds", max_peel)
)
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, output_path)
cat("Wrote RTMB retrospective output:", output_path, "\n")
