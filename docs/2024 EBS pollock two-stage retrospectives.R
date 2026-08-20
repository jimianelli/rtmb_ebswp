#!/usr/bin/env Rscript

# Nine-peel retrospective analysis that preserves the observed data lags in
# the 2024 input, fixes terminal-year fishery selectivity to the preceding
# year, and repeats the empirical-start/two-stage NonParametricPM sequence.

suppressPackageStartupMessages({
  .libPaths(c(file.path(getwd(), ".r-lib-rceattle-5.8.1"), .libPaths()))
  library(dplyr)
  library(Rceattle)
})
stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))

fit_file <- "results/canonical_pm/ebs_pollock_method_fits.rds"
scaffold_file <- "results/ebs_pollock_nonparametric_retro_9peels.rds"
output_file <- "results/canonical_pm/ebs_pollock_nonparametric_two_stage_retro_9peels.rds"
mohn_file <- "results/canonical_pm/ebs_pollock_nonparametric_two_stage_retro_9peels_mohns.csv"
peels <- 9L
n_selages_fsh <- 12L

fits <- readRDS(fit_file)
base_fit <- fits$nonparametric_pm
scaffold <- readRDS(scaffold_file)$Rceattle_list
base_endyr <- base_fit$data_list$endyr
styr <- base_fit$data_list$styr

apply_observation_lags <- function(peeled_data, base_data, endyr_peel) {
  apply_group_lag <- function(x, reference, group_var) {
    if (!nrow(x) || !nrow(reference)) return(x)
    reference_max <- reference |>
      filter(.data$Year > 0) |>
      group_by(.data[[group_var]]) |>
      summarise(max_year = max(.data$Year), .groups = "drop") |>
      mutate(
        lag = base_endyr - .data$max_year,
        cutoff = if_else(.data$lag <= peels, endyr_peel - .data$lag, .data$max_year)
      ) |>
      select(all_of(group_var), .data$cutoff)

    x |>
      left_join(reference_max, by = group_var) |>
      filter(.data$Year <= .data$cutoff | .data$Year <= 0 | is.na(.data$cutoff)) |>
      select(-.data$cutoff)
  }

  peeled_data$index_data <- apply_group_lag(
    peeled_data$index_data,
    base_data$index_data,
    "Fleet_code"
  )
  peeled_data$comp_data <- apply_group_lag(
    peeled_data$comp_data,
    base_data$comp_data,
    "Fleet_code"
  )
  peeled_data
}

empirical_start <- function(data_list) {
  yrs <- data_list$styr:data_list$endyr
  nyr <- length(yrs)
  fsh <- data_list$fleet_control$Fleet_code[
    data_list$fleet_control$Fleet_name == "Fishery"
  ]

  m0 <- fit_mod(
    data_list = data_list,
    inits = NULL,
    estimateMode = 0,
    random_rec = FALSE,
    msmMode = 0,
    initMode = "NonEquilibrium",
    M1Fun = build_M1(updateM1 = TRUE, M1_model = "fixed"),
    fit_control = fit_control(
      verbose = 0,
      phase = TRUE,
      bias_adjust_proc = 0,
      bias_adjust_obs = 0,
      comp_offset = 1e-3
    )
  )

  N <- m0$quantities$N_at_age[1, 1, , seq_len(nyr)]
  cd <- data_list$comp_data[
    data_list$comp_data$Fleet_code == fsh &
      data_list$comp_data$Year > 0 &
      data_list$comp_data$Age0_Length1 == 0,
  ]
  cc <- grep("^Comp_", colnames(cd), value = TRUE)[seq_len(data_list$nages)]
  sy <- matrix(NA_real_, nrow(cd), data_list$nages)
  for (i in seq_len(nrow(cd))) {
    yi <- which(yrs == cd$Year[i])
    if (!length(yi)) next
    pa <- as.numeric(cd[i, cc])
    pa <- pa / sum(pa, na.rm = TRUE)
    s <- pa / pmax(N[, yi], 1e-8)
    sy[i, ] <- s / max(s, na.rm = TRUE)
  }
  sel_bar <- colMeans(sy, na.rm = TRUE)[seq_len(n_selages_fsh)]
  log_sel <- log(pmax(sel_bar / max(sel_bar), 1e-3))
  log_sel <- log_sel - mean(log_sel)

  inits <- build_params(data_list)
  inits$sel_coff[1, 1, seq_len(n_selages_fsh)] <- log_sel
  inits
}

make_stage_2_map <- function(data_list, inits, endyr_peel) {
  map <- base_fit$map
  nyrs <- base_fit$data_list$endyr - styr + 1L
  nyrs_proj <- base_fit$data_list$projyr - styr + 1L
  nyrs_peel <- endyr_peel - styr + 1L

  trim_after <- list(
    rec_dev = nyrs_proj,
    log_M1_dev = nyrs_proj,
    index_q_dev = nyrs,
    log_sel_slp_dev = nyrs,
    sel_inf_dev = nyrs,
    sel_coff_dev = nyrs
  )
  for (nm in names(trim_after)) {
    x <- map$mapList[[nm]]
    year_dim <- length(dim(x))
    if (nyrs_peel < trim_after[[nm]]) {
      idx <- rep(list(TRUE), year_dim)
      idx[[year_dim]] <- (nyrs_peel + 1L):trim_after[[nm]]
      x <- do.call(`[<-`, c(list(x), idx, list(value = NA)))
    }
    map$mapList[[nm]] <- x
  }

  # The terminal fishery selectivity deviations share parameter labels with
  # the preceding year. CPUE mirrors the fishery and therefore receives the
  # same mapping. Starting values are also matched before optimization.
  terminal_i <- nyrs_peel
  previous_i <- terminal_i - 1L
  mirrored_fleets <- which(
    data_list$fleet_control$Fleet_name %in% c("Fishery", "CPUE")
  )
  for (fleet_i in mirrored_fleets) {
    map$mapList$sel_coff_dev[fleet_i, 1, , terminal_i] <-
      map$mapList$sel_coff_dev[fleet_i, 1, , previous_i]
    inits$sel_coff_dev[fleet_i, 1, , terminal_i] <-
      inits$sel_coff_dev[fleet_i, 1, , previous_i]
  }

  zero_catch <- as.matrix(
    data_list$catch_data |>
      filter(.data$Year <= base_endyr, .data$Catch == 0) |>
      mutate(Year = .data$Year - styr + 1L) |>
      select(.data$Fleet_code, .data$Year)
  )
  if (nrow(zero_catch)) {
    inits$log_F[zero_catch] <- -999
    map$mapList$log_F[zero_catch] <- NA
  }

  for (nm in names(map$mapList)) {
    map$mapFactor[[nm]] <- factor(map$mapList[[nm]])
  }
  list(map = map, inits = inits)
}

run_peel <- function(endyr_peel) {
  message("Starting two-stage peel ending in ", endyr_peel)
  peeled_data <- scaffold[[paste0("Year_", endyr_peel)]]$data_list
  peeled_data <- apply_observation_lags(
    peeled_data,
    base_fit$data_list,
    endyr_peel
  )
  peeled_data$fleet_control$Comp_distribution <- "MultinomialAFSC"

  inits <- empirical_start(peeled_data)
  stage_1_data <- peeled_data
  stage_1_data$fleet_control$Time_varying_sel <- "Off"
  ctl <- fit_control(
    verbose = 0,
    phase = TRUE,
    bias_adjust_proc = 0,
    bias_adjust_obs = 0,
    comp_offset = 1e-3
  )
  M1Fun <- build_M1(updateM1 = TRUE, M1_model = "fixed")

  stage_1 <- fit_mod(
    data_list = stage_1_data,
    inits = inits,
    file = NULL,
    estimateMode = 0,
    random_rec = FALSE,
    msmMode = 0,
    initMode = "NonEquilibrium",
    M1Fun = M1Fun,
    fit_control = ctl
  )

  stage_2_setup <- make_stage_2_map(
    peeled_data,
    stage_1$obj$env$parList(),
    endyr_peel
  )
  stage_2 <- fit_mod(
    data_list = peeled_data,
    inits = stage_2_setup$inits,
    map = stage_2_setup$map,
    file = NULL,
    estimateMode = 0,
    random_rec = FALSE,
    msmMode = 0,
    initMode = "NonEquilibrium",
    M1Fun = M1Fun,
    fit_control = ctl
  )

  list(
    terminal_year = endyr_peel,
    stage_1 = stage_1,
    stage_2 = stage_2,
    data_max_years = list(
      index = aggregate(Year ~ Fleet_code, peeled_data$index_data, max),
      composition = aggregate(Year ~ Fleet_code, peeled_data$comp_data, max)
    )
  )
}

peel_years <- (base_endyr - peels):(base_endyr - 1L)
peel_results <- lapply(peel_years, run_peel)
names(peel_results) <- paste0("Year_", peel_years)

model_list <- c(
  lapply(peel_results, `[[`, "stage_2"),
  setNames(list(base_fit), paste0("Year_", base_endyr))
)

objects <- c("biomass", "ssb", "R", "F_spp")
mohns <- bind_rows(lapply(objects, function(object) {
  relative_errors <- vapply(peel_years, function(endyr_peel) {
    i <- endyr_peel - styr + 1L
    peel <- model_list[[paste0("Year_", endyr_peel)]]$quantities[[object]][1, i]
    base <- base_fit$quantities[[object]][1, i]
    (peel - base) / base
  }, numeric(1))
  data.frame(
    Object = object,
    `Forecast year` = 0L,
    N = length(relative_errors),
    `EBS Pollock` = mean(relative_errors),
    check.names = FALSE
  )
}))

output <- list(
  Rceattle_list = model_list,
  peel_details = peel_results,
  mohns = mohns,
  specification = list(
    peels = peels,
    data_lags_preserved = TRUE,
    terminal_fishery_selectivity_equal_previous_year = TRUE,
    empirical_start_recalculated_by_peel = TRUE,
    two_stage_fit_by_peel = TRUE
    ,composition_likelihood = "MultinomialAFSC"
    ,composition_sample_sizes = "nominal, unadjusted"
    ,scaffold_role = "year-trimming scaffold only; canonical likelihood restored before every fit"
  )
)

saveRDS(output, output_file)
write.csv(mohns, mohn_file, row.names = FALSE)
message("Saved two-stage retrospectives to ", output_file)
message("Saved Mohn's rho summary to ", mohn_file)
