#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

summarize_m_profile <- function(
    profile_file = file.path(
      "results", "rceattle_5.8.1_validation", "extended_diagnostics",
      "profile_M_age3plus.rds"
    ),
    output_dir = dirname(profile_file)) {
  profile <- readRDS(profile_file)
  stopifnot(
    length(profile$Rceattle_list) == nrow(profile$grid),
    length(profile$nll) == nrow(profile$grid),
    ncol(profile$grid) == 1L
  )

  m_values <- as.numeric(profile$grid[[1]])
  min_index <- which.min(as.numeric(profile$nll))
  min_m <- m_values[min_index]

  fleet_components <- bind_rows(lapply(seq_along(profile$Rceattle_list), function(i) {
    fit <- profile$Rceattle_list[[i]]
    component_matrix <- fit$quantities$jnll_comp
    if (is.null(component_matrix) || is.null(dimnames(component_matrix))) {
      stop("Profile fit ", i, " does not contain a named jnll_comp matrix.")
    }
    as.data.frame(as.table(component_matrix), stringsAsFactors = FALSE) |>
      setNames(c("Component", "Fleet", "NLL")) |>
      mutate(
        `Age-3+ natural mortality` = m_values[i],
        .before = 1
      )
  })) |>
    group_by(Component, Fleet) |>
    mutate(
      Delta_component = NLL - NLL[`Age-3+ natural mortality` == min_m]
    ) |>
    ungroup()

  objective_check <- fleet_components |>
    group_by(`Age-3+ natural mortality`) |>
    summarize(Component_sum = sum(NLL), .groups = "drop") |>
    mutate(
      Objective = as.numeric(profile$nll),
      Difference = Component_sum - Objective
    )
  if (max(abs(objective_check$Difference)) > 1e-8) {
    stop("The jnll_comp entries do not sum to the saved profile objective.")
  }

  component_groups <- c(
    "Index data" = "Data likelihoods",
    "Catch data" = "Data likelihoods",
    "Composition data" = "Data likelihoods",
    "CAAL data" = "Data likelihoods",
    "Stomach content data" = "Data likelihoods",
    "Non-parametric selectivity" = "Selectivity penalties",
    "Selectivity deviates" = "Selectivity penalties",
    "Initial abundance deviates" = "Population penalties",
    "Recruitment deviates" = "Population penalties",
    "Stock-recruit prior" = "Population penalties",
    "Stock-recruit penalty" = "Population penalties",
    "M prior" = "Population penalties",
    "M random effects" = "Population penalties",
    "Catchability prior" = "Other penalties",
    "Catchability deviates" = "Other penalties",
    "Reference point penalties" = "Other penalties",
    "Zero n-at-age penalty" = "Other penalties",
    "Ration" = "Other penalties",
    "Ration penalties" = "Other penalties",
    "Linkage-table priors" = "Other penalties",
    "Linkage random effects" = "Other penalties"
  )

  components <- fleet_components |>
    group_by(`Age-3+ natural mortality`, Component) |>
    summarize(NLL = sum(NLL), .groups = "drop") |>
    group_by(Component) |>
    mutate(
      Delta_component = NLL - NLL[`Age-3+ natural mortality` == min_m]
    ) |>
    ungroup() |>
    mutate(Group = unname(component_groups[Component]), .before = Component)

  total <- tibble(
    `Age-3+ natural mortality` = m_values,
    Group = "Total objective",
    Component = "Total objective",
    NLL = as.numeric(profile$nll),
    Delta_component = as.numeric(profile$nll) - profile$nll[min_index]
  )

  write_csv(
    bind_rows(total, components),
    file.path(output_dir, "profile_M_age3plus_components.csv")
  )
  write_csv(
    fleet_components,
    file.path(output_dir, "profile_M_age3plus_fleet_components.csv")
  )
  write_csv(
    objective_check,
    file.path(output_dir, "profile_M_age3plus_component_check.csv")
  )

  invisible(list(
    components = bind_rows(total, components),
    fleet_components = fleet_components,
    objective_check = objective_check,
    minimum_m = min_m
  ))
}

if (sys.nframe() == 0L) {
  summarize_m_profile()
}
