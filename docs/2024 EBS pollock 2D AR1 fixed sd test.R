#!/usr/bin/env Rscript

# Test whether the 2D age-by-year AR1 fishery-selectivity sensitivity
# converges when its log standard deviation is fixed.

suppressPackageStartupMessages(library(Rceattle))

input_file <- "results/ebs_pollock_method_fits.rds"
output_file <- "results/ebs_pollock_2dar1_fixed_log_sd_test.rds"
fixed_log_sd <- -0.9702190739

fits <- readRDS(input_file)
starting_fit <- fits$ar1_2d
starting_penalized <- fits$ar1_penalized
data_list <- starting_fit$data_list

fishery_rows <- which(
  data_list$fleet_control$Fleet_name %in% c("Fishery", "CPUE")
)
# Fishery age compositions span 1964--2023. Estimate the 2D AR1 field over
# those supported years and carry the 2023 field into the 2024 terminal year.
fishery_comp_years <- sort(unique(
  data_list$comp_data$Year[
    data_list$comp_data$Fleet_code == fishery_rows[1] &
      data_list$comp_data$Age0_Length1 == 0
  ]
))
last_comp_i <- max(fishery_comp_years) - data_list$styr + 1L
terminal_i <- data_list$endyr - data_list$styr + 1L

configure_fixed_sd <- function(inits, map) {
  inits$sel_dev_log_sd[fishery_rows] <- fixed_log_sd
  map$mapList$sel_dev_log_sd[fishery_rows] <- NA
  for (fleet_i in fishery_rows) {
    inits$sel_coff_dev[fleet_i, 1, , terminal_i] <-
      inits$sel_coff_dev[fleet_i, 1, , last_comp_i]
    map$mapList$sel_coff_dev[fleet_i, 1, , terminal_i] <- NA
  }
  map$mapFactor$sel_dev_log_sd <- factor(map$mapList$sel_dev_log_sd)
  map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)
  list(inits = inits, map = map)
}

copy_matching <- function(target, source) {
  for (nm in intersect(names(target), names(source))) {
    if (identical(dim(target[[nm]]), dim(source[[nm]])) &&
        length(target[[nm]]) == length(source[[nm]])) {
      target[[nm]] <- source[[nm]]
    }
  }
  target
}

message(
  "Fitting 2D AR1 with fishery/CPUE sel_dev_log_sd fixed at ",
  fixed_log_sd,
  " (SD = ", exp(fixed_log_sd), ")."
)
message(
  "The estimated fishery/CPUE field spans ", min(fishery_comp_years),
  "--", max(fishery_comp_years), "; terminal-year values carry forward."
)

message("Stage 1: rebuilding the field as penalized effects.")
penalized_setup <- configure_fixed_sd(
  starting_penalized$estimated_params,
  starting_penalized$map
)
fixed_sd_penalized <- fit_mod(
  data_list = data_list,
  inits = penalized_setup$inits,
  map = penalized_setup$map,
  file = NULL,
  estimateMode = 0,
  random_rec = FALSE,
  random_sel = FALSE,
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

message("Stage 2: integrating the stabilized field with Laplace.")
random_inits <- copy_matching(
  starting_fit$estimated_params,
  fixed_sd_penalized$obj$env$parList()
)
random_setup <- configure_fixed_sd(random_inits, starting_fit$map)

fixed_sd_fit <- fit_mod(
  data_list = data_list,
  inits = random_setup$inits,
  map = random_setup$map,
  file = NULL,
  estimateMode = 0,
  random_rec = FALSE,
  random_sel = TRUE,
  msmMode = 0,
  initMode = "NonEquilibrium",
  M1Fun = build_M1(updateM1 = TRUE, M1_model = "fixed"),
  fit_control = fit_control(
    verbose = 1,
    phase = TRUE,
    bias_adjust_proc = 0,
    bias_adjust_obs = 0,
    comp_offset = 1e-3
  )
)

saveRDS(
  list(
    fit = fixed_sd_fit,
    penalized_fit = fixed_sd_penalized,
    fixed_log_sd = fixed_log_sd,
    fixed_sd = exp(fixed_log_sd),
    fixed_fleets = data_list$fleet_control$Fleet_name[fishery_rows]
  ),
  output_file
)

message("Saved fixed-SD test to ", output_file)
message("Status: ", fixed_sd_fit$convergence$status)
message("Maximum gradient: ", fixed_sd_fit$opt$max_gradient)
message("Convergence check: ", fixed_sd_fit$opt$Convergence_check)
