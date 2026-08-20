# Exploratory DSEM recruitment tests for the 2024 EBS pollock Rceattle model.
#
# The established assessment structure and NonParametricPM selectivity remain
# unchanged. DSEM makes recruitment deviations latent random effects and adds
# a small, prespecified set of cohort-aligned environmental paths. All candidate
# models carry the same two observed covariate series, allowing AIC comparison.

local_library <- normalizePath(".r-lib-rceattle-5.8.1-dsem", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(Rceattle)
  library(dsem)
})

stopifnot(
  packageVersion("Rceattle") == package_version("5.8.1.9000"),
  packageVersion("dsem") == package_version("3.0.0")
)

arguments <- commandArgs(trailingOnly = TRUE)
randomization_check <- "--randomized" %in% arguments
arguments <- setdiff(arguments, "--randomized")
randomization_seed <- 20260812L

baseline_directory <- "results/dsem_5.8.1"
output_directory <- if (randomization_check) {
  "results/dsem_5.8.1_randomized"
} else {
  baseline_directory
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

base_fit <- readRDS(
  file.path("results", "canonical_pm", "ebs_pollock_method_fits.rds")
)$nonparametric_pm
est <- Rceattle::read_data(
  "Data/EBS_24_pollock_m23_rceattle_full_1964-2024.xlsx"
)
est$diet_data <- NULL

# The ESP table is indexed by birth year. Both selected covariates carry the
# suffix t1 and therefore describe the following calendar year, when the cohort
# and the assessment recruitment are age 1. A zero-lag DSEM path is appropriate
# after assigning Year = cohort birth year + 1.
cohort <- read.csv(
  "/Users/jim/_mymods/pollock/ebswp_esp/data/cohort_table.csv",
  check.names = FALSE
)

standardize_observed <- function(x) {
  observed <- is.finite(x)
  answer <- rep(NA_real_, length(x))
  answer[observed] <- as.numeric(scale(x[observed]))
  answer
}

covariates <- data.frame(
  Year = cohort$year + 1L,
  ColdPool = standardize_observed(cohort$cold_pool_t1),
  AdultBiomass = standardize_observed(cohort$adult_pollock_biomass_t1)
)

if (randomization_check) {
  set.seed(randomization_seed)
  for (variable in c("ColdPool", "AdultBiomass")) {
    observed <- is.finite(covariates[[variable]])
    covariates[[variable]][observed] <- stats::rnorm(sum(observed))
  }
}

model_years <- data.frame(Year = est$styr:est$endyr)
est$env_data <- merge(model_years, covariates, by = "Year", all.x = TRUE)

write.csv(
  est$env_data,
  file.path(output_directory, "cohort_aligned_covariates.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    Rceattle_version = as.character(packageVersion("Rceattle")),
    official_base_version = "5.8.1",
    official_base_commit = "7bb694788012685495d90b806e61792e33900a50",
    dsem_patch_commit = "cb7c3f0df3b9f136085208c9711b5f8269de7352",
    dsem_version = as.character(packageVersion("dsem")),
    randomization_check = randomization_check,
    seed = if (randomization_check) randomization_seed else NA_integer_,
    distribution = if (randomization_check) "independent N(0, 1)" else "observed standardized series",
    missingness = "preserved from each observed covariate"
  ),
  file.path(output_directory, "dsem_run_design.csv"),
  row.names = FALSE
)

# Both covariates receive AR(1) latent processes and the same standardized
# observation-error scale in every model. Retaining the complete covariate
# universe makes likelihood and AIC differences attributable to recruitment
# paths rather than changing data sets.
sem_common <- c(
  "ColdPool -> ColdPool, 1, AR_ColdPool, 0",
  "AdultBiomass -> AdultBiomass, 1, AR_AdultBiomass, 0",
  "recdevs1 <-> recdevs1, 0, sigmaR1, 0.6"
)

specifications <- list(
  iid = sem_common,
  cold_pool = c(
    sem_common,
    "ColdPool -> recdevs1, 0, ColdPool_to_R, 0"
  ),
  adult_biomass = c(
    sem_common,
    "AdultBiomass -> recdevs1, 0, AdultBiomass_to_R, 0"
  ),
  combined = c(
    sem_common,
    "ColdPool -> recdevs1, 0, ColdPool_to_R, 0",
    "AdultBiomass -> recdevs1, 0, AdultBiomass_to_R, 0"
  )
)

model_display_labels <- if (randomization_check) {
  c(
    iid = "IID",
    cold_pool = "Random cold-pool substitute",
    adult_biomass = "Random SSB-proxy substitute",
    combined = "Both random substitutes"
  )
} else {
  c(
    iid = "IID",
    cold_pool = "Cold pool",
    adult_biomass = "Next-year SSB proxy",
    combined = "Cold pool + next-year SSB proxy"
  )
}

fit_one <- function(name, paths) {
  result_file <- file.path(output_directory, paste0("dsem_", name, ".rds"))
  if (file.exists(result_file)) {
    message("Using completed checkpoint: ", result_file)
    return(readRDS(result_file))
  }

  sem <- paste(paths, collapse = "\n")
  started <- Sys.time()

  # DSEM expands selected linkage parameter blocks. Begin from the branch's
  # compatible parameter structure, then transfer every baseline block whose
  # shape is unchanged. The DSEM-specific blocks are added inside fit_mod().
  compatible_inits <- Rceattle::build_params(est)
  baseline_inits <- base_fit$estimated_params
  transferable <- intersect(names(compatible_inits), names(baseline_inits))
  transferable <- transferable[vapply(
    transferable,
    function(parameter) {
      identical(dim(compatible_inits[[parameter]]), dim(baseline_inits[[parameter]])) &&
        length(compatible_inits[[parameter]]) == length(baseline_inits[[parameter]])
    },
    logical(1)
  )]
  for (parameter in transferable) {
    compatible_inits[[parameter]] <- baseline_inits[[parameter]]
  }

  model <- Rceattle::fit_mod(
    data_list = est,
    inits = compatible_inits,
    file = NULL,
    estimateMode = 0,
    random_rec = TRUE,
    msmMode = 0,
    initMode = "NonEquilibrium",
    M1Fun = Rceattle::build_M1(updateM1 = TRUE, M1_model = "fixed"),
    dsem = Rceattle::build_DSEM(
      sem = sem,
      family = dsem::gaussian_fixed_sd("identity", 0.1),
      sigmaR_prior_sd = 0.5,
      estimate_projection = FALSE
    ),
    fit_control = Rceattle::fit_control(
      verbose = 1,
      phase = TRUE,
      bias_adjust_proc = 0,
      bias_adjust_obs = 0,
      comp_offset = 1e-3
    )
  )
  attr(model, "dsem_run") <- list(
    name = name,
    sem = sem,
    started = started,
    finished = Sys.time(),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  saveRDS(model, result_file)
  model
}

requested <- arguments
if (!length(requested)) {
  requested <- names(specifications)
}
unknown <- setdiff(requested, names(specifications))
if (length(unknown)) {
  stop("Unknown DSEM configuration(s): ", paste(unknown, collapse = ", "))
}

fits <- lapply(requested, function(name) fit_one(name, specifications[[name]]))
names(fits) <- requested

# Post-processing uses the complete set of available checkpoints, even when a
# subset of configurations was requested on the command line. This permits
# interrupted or staged runs to resume without refitting completed models.
checkpoint_files <- setNames(
  file.path(output_directory, paste0("dsem_", names(specifications), ".rds")),
  names(specifications)
)
completed <- file.exists(checkpoint_files)
fits <- lapply(checkpoint_files[completed], readRDS)
names(fits) <- names(checkpoint_files)[completed]

fit_summary <- do.call(rbind, lapply(names(fits), function(name) {
  model <- fits[[name]]
  run <- attr(model, "dsem_run")
  convergence_status <- if (is.list(model$convergence)) {
    model$convergence$status
  } else {
    NA_character_
  }
  data.frame(
    model = name,
    objective = model$opt$objective,
    AIC = tryCatch(AIC(model), error = function(e) NA_real_),
    maximum_gradient = max(abs(model$obj$gr(model$opt$par))),
    positive_definite_hessian = isTRUE(model$sdrep$pdHess),
    convergence_status = convergence_status,
    elapsed_seconds = run$elapsed_seconds
  )
}))

summary_file <- file.path(output_directory, "dsem_model_summary.csv")
if (file.exists(summary_file)) {
  previous <- read.csv(summary_file)
  fit_summary <- rbind(previous[!previous$model %in% fit_summary$model, ], fit_summary)
}
fit_summary <- fit_summary[order(match(fit_summary$model, names(specifications))), ]
write.csv(fit_summary, summary_file, row.names = FALSE)
print(fit_summary, row.names = FALSE)

# Standardized path estimates and uncertainty.
path_table <- do.call(rbind, lapply(names(fits), function(name) {
  model <- fits[[name]]
  paths <- as.data.frame(summary(model$dsem))
  fixed <- as.data.frame(summary(model$sdrep, "fixed"))
  beta <- fixed[grepl("^dsem_beta_z", rownames(fixed)), , drop = FALSE]
  stopifnot(nrow(paths) == nrow(beta))
  data.frame(
    model = name,
    path = paths$path,
    lag = paths$lag,
    parameter = paths$name,
    estimate = beta$Estimate,
    standard_error = beta$`Std. Error`,
    lower_95 = beta$Estimate - 1.96 * beta$`Std. Error`,
    upper_95 = beta$Estimate + 1.96 * beta$`Std. Error`
  )
}))
write.csv(
  path_table,
  file.path(output_directory, "dsem_path_estimates.csv"),
  row.names = FALSE
)

# Quantify the change in unexplained recruitment variance relative to IID.
recruitment_scale <- path_table[path_table$parameter == "sigmaR1", ]
iid_variance <- recruitment_scale$estimate[recruitment_scale$model == "iid"]^2
recruitment_scale$unexplained_variance <- recruitment_scale$estimate^2
recruitment_scale$variance_reduction_from_iid <-
  1 - recruitment_scale$unexplained_variance / iid_variance
write.csv(
  recruitment_scale,
  file.path(output_directory, "dsem_recruitment_variance.csv"),
  row.names = FALSE
)

# Stable report figures. Figures are generated here so rendering the Quarto
# report never refits an assessment model.
suppressPackageStartupMessages(library(ggplot2))

# Covariate diagnostic: raw source series and standardized values passed to
# Rceattle. Breaks in the lines mark gaps in the annual source record.
covariate_raw <- if (randomization_check) {
  data.frame(
    Year = covariates$Year,
    `Standard-normal cold-pool substitute` = covariates$ColdPool,
    `Standard-normal SSB-proxy substitute` = covariates$AdultBiomass,
    check.names = FALSE
  )
} else {
  data.frame(
    Year = cohort$year + 1L,
    `Cold-pool extent (km2)` = cohort$cold_pool_t1,
    `Next-year SSB proxy (kt)` = cohort$adult_pollock_biomass_t1,
    check.names = FALSE
  )
}
covariate_raw <- covariate_raw |>
  tidyr::pivot_longer(
    cols = -Year,
    names_to = "Covariate",
    values_to = "Value"
  )

covariate_standardized <- est$env_data |>
  tidyr::pivot_longer(
    cols = c(ColdPool, AdultBiomass),
    names_to = "Covariate",
    values_to = "Standardized_value"
  ) |>
  dplyr::mutate(
    Covariate = dplyr::recode(
      Covariate,
      ColdPool = if (randomization_check) {
        "Standard-normal cold-pool substitute"
      } else {
        "Cold-pool extent"
      },
      AdultBiomass = if (randomization_check) {
        "Standard-normal SSB-proxy substitute"
      } else {
        "Next-year SSB proxy"
      }
    )
  )

write.csv(
  covariate_raw,
  file.path(output_directory, "dsem_covariates_raw.csv"),
  row.names = FALSE
)

raw_panel <- ggplot(covariate_raw, aes(x = Year, y = Value)) +
  geom_line(color = "#0072B2", linewidth = 0.75, na.rm = TRUE) +
  geom_point(color = "#0072B2", size = 1.5, na.rm = TRUE) +
  facet_wrap(~Covariate, ncol = 1, scales = "free_y") +
  labs(
    x = NULL,
    y = NULL,
    title = if (randomization_check) {
      "Independent standard-normal randomization series"
    } else {
      "Source series on their original scales"
    }
  ) +
  ggthemes::theme_few()

standardized_panel <- ggplot(
  covariate_standardized,
  aes(
    x = Year,
    y = Standardized_value,
    color = Covariate,
    linetype = Covariate
  )
) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_linetype_manual(values = c("longdash", "solid")) +
  labs(
    x = "Year",
    y = if (randomization_check) {
      "Standard-normal value"
    } else {
      "Standard deviations from the observed-series mean"
    },
    color = NULL,
    linetype = NULL,
    title = "Standardized covariates supplied to DSEM"
  ) +
  ggthemes::theme_few() +
  theme(legend.position = "bottom")

covariate_figure <- patchwork::wrap_plots(
  raw_panel,
  standardized_panel,
  ncol = 1,
  heights = c(1.4, 1)
) +
  patchwork::plot_annotation(
    caption = if (randomization_check) {
      paste0(
        "Independent N(0, 1) values retain the observed missing-year patterns; ",
        "seed = ", randomization_seed, "."
      )
    } else {
      paste0(
        "Cold-pool observations: ", sum(is.finite(est$env_data$ColdPool)),
        "; next-year SSB-proxy observations: ",
        sum(is.finite(est$env_data$AdultBiomass)), "."
      )
    }
  )

ggsave(
  file.path(output_directory, "dsem_covariates.png"),
  covariate_figure,
  width = 9,
  height = 9,
  dpi = 180
)

fit_summary$delta_AIC <- fit_summary$AIC - min(fit_summary$AIC)
fit_summary$label <- unname(model_display_labels[fit_summary$model])

aic_figure <- ggplot(fit_summary, aes(x = reorder(label, AIC), y = delta_AIC)) +
  geom_col(fill = "#176b87", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", delta_AIC)), hjust = -0.15) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression(Delta * "AIC")) +
  ggthemes::theme_few()
ggsave(
  file.path(output_directory, "dsem_aic_comparison.png"),
  aic_figure,
  width = 7,
  height = 4.5,
  dpi = 180
)

years <- est$styr:est$endyr
recruitment <- do.call(rbind, lapply(names(fits), function(name) {
  values <- as.numeric(fits[[name]]$quantities$R[1, seq_along(years)])
  data.frame(Year = years, Recruitment = values, model = name)
}))
recruitment$model <- factor(
  recruitment$model,
  levels = names(specifications),
  labels = unname(model_display_labels[names(specifications)])
)
write.csv(
  recruitment,
  file.path(output_directory, "dsem_recruitment_trajectories.csv"),
  row.names = FALSE
)

recruitment_figure <- ggplot(
  recruitment,
  aes(x = Year, y = Recruitment, color = model, linetype = model)
) +
  geom_line(linewidth = 0.85, alpha = 0.95) +
  scale_color_manual(values = c("#000000", "#0072B2", "#D55E00", "#009E73")) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotdash")) +
  labs(
    x = "Year",
    y = "Age-1 recruitment",
    color = "Model",
    linetype = "Model"
  ) +
  ggthemes::theme_few() +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(1.8, "cm")
  )
ggsave(
  file.path(output_directory, "dsem_recruitment_comparison.png"),
  recruitment_figure,
  width = 9,
  height = 5,
  dpi = 180
)

# Directed-edge coefficient plot. Variance paths use the bidirectional `<->`
# notation and are omitted because they represent scales rather than effects
# from one variable to another.
edge_coefficients <- path_table |>
  dplyr::filter(grepl(" -> ", path, fixed = TRUE)) |>
  dplyr::mutate(
    model_label = factor(
      model,
      levels = names(specifications),
      labels = unname(model_display_labels[names(specifications)])
    ),
    edge = dplyr::recode(
      parameter,
      AR_ColdPool = if (randomization_check) {
        "Random cold-pool substitute(t-1) -> value(t)"
      } else {
        "Cold pool(t-1) -> cold pool(t)"
      },
      AR_AdultBiomass = if (randomization_check) {
        "Random SSB-proxy substitute(t-1) -> value(t)"
      } else {
        "Next-year SSB proxy(t-1) -> proxy(t)"
      },
      ColdPool_to_R = if (randomization_check) {
        "Random cold-pool substitute -> age-1 recruitment"
      } else {
        "Cold pool -> age-1 recruitment"
      },
      AdultBiomass_to_R = if (randomization_check) {
        "Random SSB-proxy substitute -> age-1 recruitment"
      } else {
        "Next-year SSB proxy -> age-1 recruitment"
      }
    ),
    edge_type = ifelse(
      grepl("_to_R$", parameter),
      "Recruitment effect",
      "Environmental persistence"
    )
  )

write.csv(
  edge_coefficients,
  file.path(output_directory, "dsem_edge_coefficients.csv"),
  row.names = FALSE
)

edge_figure <- ggplot(
  edge_coefficients,
  aes(x = estimate, y = edge, color = edge_type)
) +
  geom_vline(xintercept = 0, color = "grey55", linetype = "dashed") +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0,
    linewidth = 0.65
  ) +
  geom_point(size = 2.4) +
  facet_wrap(~model_label, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Environmental persistence" = "#176b87",
      "Recruitment effect" = "#b44b27"
    )
  ) +
  labs(
    x = "Standardized coefficient (95% interval)",
    y = NULL,
    color = "Edge type"
  ) +
  ggthemes::theme_few() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(output_directory, "dsem_edge_coefficients.png"),
  edge_figure,
  width = 10,
  height = 7,
  dpi = 180
)

capture.output(
  sessionInfo(),
  file = file.path(output_directory, "session_info.txt")
)
