#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

rtmb_root <- normalizePath(getwd(), mustWork = TRUE)
afscosa_root <- Sys.getenv(
  "AFSCOSA_ROOT",
  "/Users/jim/_mymods/noaa-afsc/afscOSA"
)
afscosa_root <- normalizePath(afscosa_root, mustWork = TRUE)
devtools::load_all(afscosa_root, quiet = TRUE)

out_dir <- file.path(rtmb_root, "analysis", "output", "osa_selectivity_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_paths <- c(
  "Base coefficients" = file.path(
    rtmb_root, "analysis", "output", "fishery_sel_forms", "fishery_sel_form_0.rds"
  ),
  "Time-age varying double logistic (20% CV)" = file.path(
    rtmb_root, "analysis", "output", "double_logistic_experiments",
    "cv0p2_no_old_age_cap_stage2_random.rds"
  ),
  "2D AR1" = file.path(
    rtmb_root, "analysis", "output", "fishery_sel_forms", "fishery_sel_form_5.rds"
  )
)

missing_runs <- run_paths[!file.exists(run_paths)]
if (length(missing_runs)) {
  stop("Missing selectivity run(s): ", paste(missing_runs, collapse = ", "))
}
runs <- lapply(run_paths, readRDS)

gear_specs <- list(
  fishery = list(
    title = "Fishery age compositions",
    obs = "oac_fsh", pred = "phat_fsh", n = "sam_fsh", years = "yrs_fsh_data",
    bins = 1:15
  ),
  bts = list(
    title = "Bottom trawl survey age compositions",
    obs = "oac_bts", pred = "phat_bts", n = "sam_bts", years = "yrs_bts_data",
    bins = 1:15
  ),
  ats = list(
    title = "Acoustic-trawl survey age compositions",
    obs = "oac_ats", pred = "phat_ats", n = "sam_ats", years = "yrs_ats_data",
    bins = 2:15
  )
)

extract_report <- function(x) x$report

validate_common_inputs <- function(reports, spec, gear) {
  reference <- reports[[1]]
  fields <- c(spec$obs, spec$n, spec$years)
  for (field in fields) {
    for (i in seq_along(reports)[-1]) {
      # Compact random-effects checkpoints do not duplicate common inputs.
      if (is.null(reports[[i]][[field]])) next
      same <- isTRUE(all.equal(reference[[field]], reports[[i]][[field]],
                               tolerance = 1e-10))
      if (!same) stop("Input mismatch for ", gear, " field ", field)
    }
  }
}

osa_summary_row <- function(x, gear, approach) {
  tibble(
    gear = gear,
    approach = approach,
    n_residuals = nrow(x$res),
    sdnr = sd(x$res$resid),
    mean = mean(x$res$resid),
    lower_2.5 = unname(quantile(x$res$resid, 0.025)),
    upper_97.5 = unname(quantile(x$res$resid, 0.975)),
    proportion_abs_gt_2 = mean(abs(x$res$resid) > 2),
    proportion_abs_gt_3 = mean(abs(x$res$resid) > 3),
    aggregate_ess = x$agg$ESS[1],
    aggregate_iss = x$agg$ISS[1]
  )
}

reports <- lapply(runs, extract_report)
all_outputs <- list()
all_summaries <- list()

for (gear in names(gear_specs)) {
  spec <- gear_specs[[gear]]
  validate_common_inputs(reports, spec, gear)
  common_report <- reports[[1]]
  gear_outputs <- imap(reports, function(report, approach) {
    bins <- spec$bins
    obs <- common_report[[spec$obs]][, bins, drop = FALSE]
    expected <- report[[spec$pred]][, bins, drop = FALSE]
    sample_size <- common_report[[spec$n]]
    years <- common_report[[spec$years]]
    rounded_n <- rowSums(round(sample_size * obs / rowSums(obs), 0))
    keep <- is.finite(rounded_n) & rounded_n >= 1
    if (any(!keep)) {
      message(
        gear, ": omitting composition year(s) with zero rounded counts: ",
        paste(years[!keep], collapse = ", ")
      )
    }
    afscOSA::run_osa(
      obs = obs[keep, , drop = FALSE],
      exp = expected[keep, , drop = FALSE],
      N = sample_size[keep],
      fleet = approach,
      index = bins,
      years = years[keep],
      index_label = "Age",
      seed = 202508L,
      nonfinite_action = "truncate",
      nonfinite_limit = 6
    )
  })
  all_outputs[[gear]] <- gear_outputs
  all_summaries[[gear]] <- imap_dfr(
    gear_outputs,
    function(x, approach) osa_summary_row(x, gear, approach)
  )

  grDevices::pdf(file = NULL)
  figure <- afscOSA::plot_osa(
    gear_outputs,
    plot = TRUE,
    add_agg_CI = TRUE,
    add_sdnr_CI = TRUE,
    add_QQ_quantiles = TRUE,
    use_agg_proportions = TRUE,
    figheight = 13,
    figwidth = 16
  )
  grDevices::dev.off()
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("osa_", gear, "_selectivity_comparison.png")),
    plot = figure, width = 16, height = 13, units = "in", dpi = 220,
    bg = "white"
  )
}

summary <- bind_rows(all_summaries)
write_csv(summary, file.path(out_dir, "osa_selectivity_summary.csv"))
saveRDS(
  list(
    afscOSA_commit = system2(
      "git", c("-C", afscosa_root, "rev-parse", "HEAD"), stdout = TRUE
    ),
    seed = 202508L,
    run_paths = run_paths,
    run_md5 = unname(tools::md5sum(run_paths)),
    form2_configuration = list(
      hierarchical_stage = runs[[2]]$hierarchical_stage,
      process_cv = runs[[2]]$cv,
      old_age_cap = runs[[2]]$old_age_cap,
      random_effects = runs[[2]]$random_effects
    ),
    outputs = all_outputs,
    summary = summary
  ),
  file.path(out_dir, "osa_selectivity_outputs.rds")
)
print(summary, n = Inf)
