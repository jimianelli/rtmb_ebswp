#!/usr/bin/env Rscript

# Diagnose the 2012--2014 cohort revisions across the 2020-to-2021
# Rceattle retrospective transition. All controlled refits start from the
# converged 2021-peel solution and retain its parameter map. They isolate local
# data and selectivity effects within the same optimization mode; they are
# diagnostics rather than alternative assessment configurations.

suppressPackageStartupMessages({
  .libPaths(c(file.path(getwd(), ".r-lib-rceattle-5.8.1"), .libPaths()))
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(readr)
  library(Rceattle)
  library(tidyr)
})

stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))

retro_file <- file.path(
  "results", "canonical_pm",
  "ebs_pollock_nonparametric_two_stage_retro_9peels.rds"
)
output_dir <- file.path(
  "results", "canonical_pm", "cohort_revision_diagnostics"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

retro <- readRDS(retro_file)
models <- retro$Rceattle_list
model_2020 <- models$Year_2020
model_2021 <- models$Year_2021
cohorts <- 2012:2014

ctl <- fit_control(
  verbose = 0,
  phase = TRUE,
  bias_adjust_proc = 0,
  bias_adjust_obs = 0,
  comp_offset = 1e-3
)
M1_fun <- build_M1(updateM1 = TRUE, M1_model = "fixed")

fit_variant <- function(data_list, inits = model_2021$obj$env$parList(),
                        map = model_2021$map) {
  suppressWarnings(fit_mod(
    data_list = data_list,
    inits = inits,
    map = map,
    file = NULL,
    estimateMode = 0,
    random_rec = FALSE,
    msmMode = 0,
    initMode = "NonEquilibrium",
    M1Fun = M1_fun,
    fit_control = ctl
  ))
}

extract_state <- function(name, fit, experiment, reference = model_2021) {
  years <- as.integer(colnames(fit$quantities$R))
  reference_years <- as.integer(colnames(reference$quantities$R))
  values <- tibble(
    Experiment = experiment,
    Scenario = name,
    Maximum_gradient = max(abs(fit$obj$gr(fit$opt$par))),
    Objective = fit$opt$objective,
    R2012 = fit$quantities$R[1, match(2012, years)],
    R2013 = fit$quantities$R[1, match(2013, years)],
    R2014 = fit$quantities$R[1, match(2014, years)],
    SSB2015 = fit$quantities$ssb[1, match(2015, years)],
    SSB2018 = fit$quantities$ssb[1, match(2018, years)],
    SSB2020 = fit$quantities$ssb[1, match(2020, years)]
  )
  for (quantity in c("R2012", "R2013", "R2014", "SSB2015", "SSB2018", "SSB2020")) {
    reference_value <- switch(
      quantity,
      R2012 = reference$quantities$R[1, match(2012, reference_years)],
      R2013 = reference$quantities$R[1, match(2013, reference_years)],
      R2014 = reference$quantities$R[1, match(2014, reference_years)],
      SSB2015 = reference$quantities$ssb[1, match(2015, reference_years)],
      SSB2018 = reference$quantities$ssb[1, match(2018, reference_years)],
      SSB2020 = reference$quantities$ssb[1, match(2020, reference_years)]
    )
    values[[paste0(quantity, "_percent_vs_all_2021")]] <-
      100 * (values[[quantity]] / reference_value - 1)
  }
  values
}

# Cohort estimates across all peels.
cohort_revisions <- bind_rows(lapply(names(models), function(model_name) {
  fit <- models[[model_name]]
  years <- as.integer(colnames(fit$quantities$R))
  tibble(
    Terminal_year = as.integer(sub("Year_", "", model_name)),
    Cohort = cohorts,
    Age1_abundance = as.numeric(fit$quantities$R[1, match(cohorts, years)])
  )
}))
write_csv(cohort_revisions, file.path(output_dir, "cohort_revisions.csv"))

# Approximate 95% uncertainty ribbons, matching the historical ADMB display
# convention of estimate +/- 2 standard deviations.
ssb_uncertainty <- bind_rows(lapply(names(models), function(model_name) {
  fit <- models[[model_name]]
  terminal_year <- as.integer(sub("Year_", "", model_name))
  years <- fit$data_list$styr:terminal_year
  report_summary <- summary(fit$sdrep, "report")
  ssb_rows <- which(rownames(report_summary) == "ssb")[seq_along(years)]
  tibble(
    Terminal_year = terminal_year,
    Year = years,
    Estimate = as.numeric(fit$quantities$ssb[1, seq_along(years)]),
    Standard_error = report_summary[ssb_rows, "Std. Error"],
    Lower_2SD = pmax(0, Estimate - 2 * Standard_error),
    Upper_2SD = Estimate + 2 * Standard_error
  )
}))
write_csv(ssb_uncertainty, file.path(output_dir, "ssb_retrospective_uncertainty.csv"))

# Inventory observations entering between the 2020 and 2021 peels.
new_rows <- function(old, new, source) {
  key_old <- paste(old$Fleet_name, old$Year)
  added <- new |>
    filter(!paste(.data$Fleet_name, .data$Year) %in% key_old)

  if (source == "Index") {
    return(added |>
      transmute(
        Source = source,
        Fleet = .data$Fleet_name,
        Year = .data$Year,
        Value = .data$Observation,
        Value_definition = "Observation"
      ))
  }

  added |>
    transmute(
      Source = source,
      Fleet = .data$Fleet_name,
      Year = .data$Year,
      Value = .data$Sample_size,
      Value_definition = "Sample size"
    )
}
new_data_inventory <- bind_rows(
  new_rows(
    model_2020$data_list$index_data,
    model_2021$data_list$index_data,
    "Index"
  ),
  new_rows(
    model_2020$data_list$comp_data,
    model_2021$data_list$comp_data,
    "Composition"
  )
)
write_csv(new_data_inventory, file.path(output_dir, "new_data_inventory.csv"))

remove_block <- list(
  no_BTS_index_2021 = function(data) {
    data$index_data <- filter(
      data$index_data,
      !(.data$Fleet_name == "BTS" & .data$Year == 2021)
    )
    data
  },
  no_BTS_comp_2021 = function(data) {
    data$comp_data <- filter(
      data$comp_data,
      !(.data$Fleet_name == "BTS" & .data$Year == 2021)
    )
    data
  },
  no_AVO_index_2021 = function(data) {
    data$index_data <- filter(
      data$index_data,
      !(.data$Fleet_name == "AVO" & .data$Year == 2021)
    )
    data
  },
  no_Fishery_comp_2020 = function(data) {
    data$comp_data <- filter(
      data$comp_data,
      !(.data$Fleet_name == "Fishery" & .data$Year == 2020)
    )
    data
  },
  no_BTS1_index_2021 = function(data) {
    data$index_data <- filter(
      data$index_data,
      !(.data$Fleet_name == "BTS_1" & .data$Year == 2021)
    )
    data
  }
)

block_deletion <- list(extract_state(
  "All 2021-peel data", model_2021, "Block deletion"
))
for (scenario in names(remove_block)) {
  message("Block deletion: ", scenario)
  fit <- fit_variant(remove_block[[scenario]](model_2021$data_list))
  block_deletion[[length(block_deletion) + 1L]] <- extract_state(
    scenario, fit, "Block deletion"
  )
}
block_deletion <- bind_rows(block_deletion)
write_csv(block_deletion, file.path(output_dir, "block_deletion_results.csv"))

# Cumulative additions use the 2021 state/process year but begin with the 2020
# observation availability. The ordering is explicit because nonlinear joint
# effects need not be additive.
full_data <- model_2021$data_list
base_availability <- full_data
for (scenario in c(
  "no_BTS_index_2021", "no_BTS_comp_2021", "no_AVO_index_2021",
  "no_Fishery_comp_2020"
)) {
  base_availability <- remove_block[[scenario]](base_availability)
}

add_rows <- function(data, source_data, source, fleet, year) {
  if (source == "index") {
    data$index_data <- bind_rows(
      data$index_data,
      filter(source_data$index_data, .data$Fleet_name == fleet, .data$Year == year)
    ) |>
      arrange(.data$Fleet_code, .data$Year)
  } else {
    data$comp_data <- bind_rows(
      data$comp_data,
      filter(source_data$comp_data, .data$Fleet_name == fleet, .data$Year == year)
    ) |>
      arrange(.data$Fleet_code, .data$Year)
  }
  data
}

cumulative_steps <- list(
  "2020 data availability" = base_availability
)
working <- base_availability
working <- add_rows(working, full_data, "composition", "Fishery", 2020)
cumulative_steps[["Add 2020 fishery composition"]] <- working
working <- add_rows(working, full_data, "index", "BTS", 2021)
cumulative_steps[["Add 2021 BTS index"]] <- working
working <- add_rows(working, full_data, "composition", "BTS", 2021)
cumulative_steps[["Add 2021 BTS composition"]] <- working
working <- add_rows(working, full_data, "index", "AVO", 2021)
cumulative_steps[["Add 2021 AVO index"]] <- working

cumulative_results <- list()
for (scenario in names(cumulative_steps)) {
  message("Cumulative addition: ", scenario)
  fit <- fit_variant(cumulative_steps[[scenario]])
  cumulative_results[[length(cumulative_results) + 1L]] <- extract_state(
    scenario, fit, "Cumulative addition"
  )
}
cumulative_results <- bind_rows(cumulative_results) |>
  mutate(Step = row_number())
write_csv(cumulative_results, file.path(output_dir, "cumulative_addition_results.csv"))

# Hold the realized fishery/CPUE or BTS selectivity parameter fields at their
# 2020-peel values, then refit all other parameters to the 2021 data.
fit_fixed_selectivity <- function(fleet_group) {
  inits <- model_2021$obj$env$parList()
  inits_2020 <- model_2020$obj$env$parList()
  map <- model_2021$map

  if (fleet_group == "Fishery and CPUE") {
    fleets <- c(1L, 7L)
    inits$sel_coff[fleets, , ] <- inits_2020$sel_coff[fleets, , ]
    map$mapList$sel_coff[fleets, , ] <- NA
    map$mapFactor$sel_coff <- factor(map$mapList$sel_coff)
    inits$sel_coff_dev[fleets, , , ] <- inits_2020$sel_coff_dev[fleets, , , ]
    map$mapList$sel_coff_dev[fleets, , , ] <- NA
    map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)
  } else if (fleet_group == "BTS") {
    fleet <- 3L
    for (parameter in c("log_sel_slp", "sel_inf")) {
      inits[[parameter]][, fleet, ] <- inits_2020[[parameter]][, fleet, ]
      map$mapList[[parameter]][, fleet, ] <- NA
      map$mapFactor[[parameter]] <- factor(map$mapList[[parameter]])
    }
    for (parameter in c("log_sel_slp_dev", "sel_inf_dev")) {
      inits[[parameter]][, fleet, , ] <- inits_2020[[parameter]][, fleet, , ]
      map$mapList[[parameter]][, fleet, , ] <- NA
      map$mapFactor[[parameter]] <- factor(map$mapList[[parameter]])
    }
  }
  fit_variant(model_2021$data_list, inits = inits, map = map)
}

selectivity_controls <- list(extract_state(
  "All selectivity re-estimated", model_2021, "Selectivity control"
))
for (fleet_group in c("Fishery and CPUE", "BTS")) {
  message("Fixed selectivity: ", fleet_group)
  fit <- fit_fixed_selectivity(fleet_group)
  selectivity_controls[[length(selectivity_controls) + 1L]] <- extract_state(
    paste0(fleet_group, " fixed at 2020 peel"),
    fit,
    "Selectivity control"
  )
}
selectivity_controls <- bind_rows(selectivity_controls)
write_csv(
  selectivity_controls,
  file.path(output_dir, "selectivity_control_results.csv")
)

# Observed-versus-predicted cohort traces for the 2021 peel.
comp_data <- model_2021$data_list$comp_data
cohort_composition <- bind_rows(lapply(seq_len(nrow(comp_data)), function(row) {
  bind_rows(lapply(cohorts, function(cohort) {
    age <- comp_data$Year[row] - cohort + 1L
    if (age < 1L || age > 15L) return(NULL)
    tibble(
      Fleet = comp_data$Fleet_name[row],
      Year = comp_data$Year[row],
      Cohort = cohort,
      Age = age,
      Sample_size = comp_data$Sample_size[row],
      Observed = model_2021$quantities$comp_obs[row, age],
      Predicted = model_2021$quantities$comp_hat[row, age],
      Residual = Observed - Predicted,
      Weighted_absolute_difference = Sample_size * abs(Residual)
    )
  }))
}))
write_csv(
  cohort_composition,
  file.path(output_dir, "cohort_composition_trace.csv")
)

# Recent index fits and 2020-to-2021 realized selectivity changes.
index_data <- model_2021$data_list$index_data
recent_index_fit <- index_data |>
  mutate(
    Predicted = as.numeric(model_2021$quantities$index_hat),
    Observed_to_predicted = .data$Observation / .data$Predicted
  ) |>
  filter(.data$Fleet_name %in% c("BTS", "AVO"), .data$Year >= 2017) |>
  select(.data$Fleet_name, .data$Year, .data$Observation, .data$Predicted,
         .data$Observed_to_predicted)
write_csv(recent_index_fit, file.path(output_dir, "recent_index_fit.csv"))

selectivity_change <- bind_rows(lapply(c("Fishery", "BTS", "ATS"), function(fleet) {
  year_index_2020 <- match(
    "2020", dimnames(model_2020$quantities$sel_at_age)[[4]]
  )
  year_index_2021 <- match(
    "2020", dimnames(model_2021$quantities$sel_at_age)[[4]]
  )
  old <- model_2020$quantities$sel_at_age[fleet, 1, , year_index_2020]
  new <- model_2021$quantities$sel_at_age[fleet, 1, , year_index_2021]
  tibble(
    Fleet = fleet,
    Age = seq_along(old),
    Selectivity_2020_peel = as.numeric(old),
    Selectivity_2021_peel = as.numeric(new),
    Ratio_2021_to_2020 = as.numeric(new / old)
  )
}))
write_csv(selectivity_change, file.path(output_dir, "selectivity_change.csv"))

# Likelihood-component changes are descriptive because the 2021 peel includes
# additional observations. They identify where the refitted objective moved but
# are not a like-for-like likelihood-ratio test.
likelihood_components <- bind_rows(lapply(
  c("Year_2020", "Year_2021"),
  function(model_name) {
    fit <- models[[model_name]]
    as.data.frame(as.table(fit$quantities$jnll_comp)) |>
      rename(Component = .data$Var1, Fleet = .data$Var2, NLL = .data$Freq) |>
      mutate(Terminal_year = as.integer(sub("Year_", "", model_name)))
  }
)) |>
  filter(is.finite(.data$NLL), abs(.data$NLL) > 1e-10)
write_csv(
  likelihood_components,
  file.path(output_dir, "likelihood_components_2020_2021.csv")
)

# Figures ---------------------------------------------------------------------
theme_set(ggthemes::theme_few())

p_cohort <- ggplot(
  cohort_revisions,
  aes(.data$Terminal_year, .data$Age1_abundance, color = factor(.data$Cohort))
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Retrospective terminal year",
    y = "Estimated age-1 abundance",
    color = "Cohort"
  )
ggsave(
  file.path(output_dir, "cohort_revisions.png"),
  p_cohort, width = 8, height = 5, dpi = 200
)

p_uncertainty <- ggplot(
  filter(ssb_uncertainty, .data$Year >= 2000),
  aes(.data$Year, .data$Estimate, group = factor(.data$Terminal_year),
      color = factor(.data$Terminal_year), fill = factor(.data$Terminal_year))
) +
  geom_ribbon(
    aes(ymin = .data$Lower_2SD, ymax = .data$Upper_2SD),
    alpha = 0.08, color = NA
  ) +
  geom_line(linewidth = 0.7) +
  scale_color_viridis_d(option = "C", end = 0.9, direction = -1) +
  scale_fill_viridis_d(option = "C", end = 0.9, direction = -1) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = NULL,
    y = "Spawning biomass (kt)",
    color = "Terminal year",
    fill = "Terminal year"
  ) +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "ssb_retrospective_uncertainty.png"),
  p_uncertainty, width = 10, height = 7, dpi = 200
)

block_plot <- block_deletion |>
  filter(.data$Scenario != "All 2021-peel data") |>
  select(
    .data$Scenario,
    R2013 = .data$R2013_percent_vs_all_2021,
    R2014 = .data$R2014_percent_vs_all_2021,
    SSB2018 = .data$SSB2018_percent_vs_all_2021,
    SSB2020 = .data$SSB2020_percent_vs_all_2021
  ) |>
  pivot_longer(-.data$Scenario, names_to = "Quantity", values_to = "Percent") |>
  mutate(
    Scenario = recode(
      .data$Scenario,
      no_BTS_index_2021 = "Remove 2021 BTS index",
      no_BTS_comp_2021 = "Remove 2021 BTS composition",
      no_AVO_index_2021 = "Remove 2021 AVO index",
      no_Fishery_comp_2020 = "Remove 2020 fishery composition",
      no_BTS1_index_2021 = "Remove 2021 BTS age-1 index"
    )
  )
p_blocks <- ggplot(
  block_plot,
  aes(.data$Percent, .data$Scenario, fill = .data$Quantity)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  labs(
    x = "Percent change from the all-data 2021 peel",
    y = NULL,
    fill = "Estimate"
  ) +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "block_deletion_effects.png"),
  p_blocks, width = 10, height = 5.5, dpi = 200
)

composition_plot_data <- cohort_composition |>
  filter(
    .data$Fleet %in% c("Fishery", "BTS", "ATS"),
    .data$Year >= 2013,
    .data$Year <= 2021
  ) |>
  pivot_longer(
    c(.data$Observed, .data$Predicted),
    names_to = "Series", values_to = "Proportion"
  )
p_composition <- ggplot(
  composition_plot_data,
  aes(.data$Year, .data$Proportion, color = .data$Series,
      linetype = .data$Series)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_grid(rows = vars(.data$Cohort), cols = vars(.data$Fleet), scales = "free_y") +
  scale_color_manual(values = c(Observed = "black", Predicted = "#D55E00")) +
  scale_linetype_manual(values = c(Observed = "solid", Predicted = "dashed")) +
  labs(
    x = "Observation year",
    y = "Cohort proportion in age composition",
    color = NULL,
    linetype = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "cohort_composition_trace.png"),
  p_composition, width = 11, height = 8, dpi = 200
)

write_lines(
  c(
    paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Rceattle version:", as.character(packageVersion("Rceattle"))),
    paste("Source retrospective:", retro_file),
    "Controlled refits start from the converged 2021-peel solution.",
    "Cumulative-addition order: 2020 fishery composition, 2021 BTS index, 2021 BTS composition, 2021 AVO index."
  ),
  file.path(output_dir, "lineage.txt")
)

message("Wrote cohort-revision diagnostics to ", output_dir)
