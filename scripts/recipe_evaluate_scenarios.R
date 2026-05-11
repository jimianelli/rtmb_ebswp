#!/usr/bin/env Rscript

# Recipe: Evaluate fishery selectivity scenarios (forms 0, 2, 5)
#
# What this script does:
#   1) Ensures POLLOCK_ROOT is set and looks valid
#   2) Runs the scenario fits via R/run_fishery_selectivity_forms.R
#   3) Loads the resulting *.rds files and performs quick diagnostics:
#        - convergence/max gradient
#        - objective deltas
#        - check that age comps (oac/phat) sum to 1 within each year (by fleet)
#        - catch observed vs predicted summary
#   4) Renders the comparison report to reporting/ and copies to docs/
#   5) Prints a git commit/push checklist (does NOT push automatically)
#
# Usage:
#   POLLOCK_ROOT=/path/to/pollock Rscript scripts/recipe_evaluate_scenarios.R
#
# Notes:
#   - This assumes you have already set up rtmb_ebswp dependencies.
#   - This script does not modify model code; it only runs and validates outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

repo_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "R", "run_fishery_selectivity_forms.R"))) {
  # If launched from scripts/, climb up
  repo_root2 <- normalizePath(file.path(repo_root, ".."), mustWork = TRUE)
  if (file.exists(file.path(repo_root2, "R", "run_fishery_selectivity_forms.R"))) {
    repo_root <- repo_root2
  }
}

cat("Repo root:", repo_root, "\n")
setwd(repo_root)

pollock_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
if (is.na(pollock_root) || pollock_root == "") {
  stop("POLLOCK_ROOT is not set. Example:\n  POLLOCK_ROOT=/Users/you/workspace/pollock Rscript scripts/recipe_evaluate_scenarios.R")
}

pollock_root <- normalizePath(pollock_root, mustWork = TRUE)
cat("POLLOCK_ROOT:", pollock_root, "\n")

# minimal structural check
need <- c(
  file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.par"),
  file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.rep"),
  file.path(pollock_root, "admb", "runs", "for_rtmb", "pm.tpl")
)
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop("POLLOCK_ROOT does not look like the expected pollock workspace. Missing:\n", paste0("- ", missing, collapse = "\n"))
}

# ---- 1) Run scenarios ----
cat("\n== Running scenario fits ==\n")
cmd <- sprintf("POLLOCK_ROOT=%s Rscript R/run_fishery_selectivity_forms.R", shQuote(pollock_root))
status <- system(cmd)
if (!identical(status, 0L)) stop("Scenario fit script failed with status=", status)

# ---- 2) Load outputs and summarize ----
results_dir <- file.path("analysis", "output", "fishery_sel_forms")
forms <- c(0, 2, 5)
paths <- file.path(results_dir, sprintf("fishery_sel_form_%d.rds", forms))
if (any(!file.exists(paths))) {
  stop("Missing one or more scenario output files:\n", paste0("- ", paths[!file.exists(paths)], collapse = "\n"))
}

runs <- lapply(seq_along(forms), function(i) {
  x <- readRDS(paths[i])
  list(
    form = forms[i],
    path = paths[i],
    label = x$label %||% paste0("form ", forms[i]),
    objective = x$objective,
    convergence = x$convergence,
    max_gradient = x$max_gradient,
    seconds = x$seconds,
    report = x$report
  )
})

# helper
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

summary_tbl <- tibble(
  fishery_sel_form = vapply(runs, `[[`, integer(1), "form"),
  label = vapply(runs, `[[`, character(1), "label"),
  objective = vapply(runs, `[[`, numeric(1), "objective"),
  convergence = vapply(runs, `[[`, integer(1), "convergence"),
  max_gradient = vapply(runs, `[[`, numeric(1), "max_gradient"),
  seconds = vapply(runs, `[[`, numeric(1), "seconds")
) |>
  arrange(fishery_sel_form) |>
  mutate(delta_objective = objective - objective[fishery_sel_form == 0][1])

cat("\n== Optimization summary ==\n")
print(summary_tbl)

# ---- 3) Diagnostics: age-comp sums ----
cat("\n== Age-comp sum-to-1 checks (within year) ==\n")

check_comp <- function(mat, name) {
  if (is.null(mat)) return(tibble(component = name, ok = NA, min = NA, max = NA))
  rs <- rowSums(mat)
  tibble(
    component = name,
    ok = all(is.finite(rs)) && max(abs(rs - 1)) < 1e-6,
    min = min(rs),
    max = max(rs)
  )
}

for (r in runs) {
  rep <- r$report
  cat("\nScenario form=", r$form, " (", r$label, ")\n", sep = "")

  # Observed
  print(bind_rows(
    check_comp(rep$oac_fsh, "oac_fsh"),
    check_comp(rep$oac_bts, "oac_bts"),
    check_comp(rep$oac_ats, "oac_ats")
  ))

  # Expected
  print(bind_rows(
    check_comp(rep$phat_fsh, "phat_fsh"),
    check_comp(rep$phat_bts, "phat_bts"),
    check_comp(rep$phat_ats, "phat_ats")
  ))
}

# ---- 4) Catch fit quick check ----
cat("\n== Catch fit quick check ==\n")
for (r in runs) {
  rep <- r$report
  obs <- rep$obs_catch
  pred <- rep$pred_catch
  if (is.null(obs) || is.null(pred)) {
    cat("form=", r$form, ": obs_catch/pred_catch not found in report\n", sep = "")
    next
  }
  yrs <- as.integer(names(obs) %||% names(pred) %||% seq_along(obs))
  df <- tibble(year = yrs, obs = as.numeric(obs), pred = as.numeric(pred))
  df <- df |> filter(is.finite(obs), is.finite(pred))
  cat("form=", r$form, ": corr(obs,pred)=", round(cor(df$obs, df$pred), 3),
      ", mean(|log(obs/pred)|)=", round(mean(abs(log((df$obs + 1e-12)/(df$pred + 1e-12)))), 3), "\n", sep = "")
}

# ---- 5) Render report ----
cat("\n== Render comparison report ==\n")
qmd <- file.path("reporting", "fishery_sel_form_comparison.qmd")
if (!file.exists(qmd)) stop("Missing: ", qmd)

status <- system(sprintf("quarto render %s --no-cache", shQuote(qmd)))
if (!identical(status, 0L)) stop("quarto render failed")

# copy to docs for Pages
html_out <- file.path("reporting", "fishery_sel_form_comparison.html")
if (!file.exists(html_out)) stop("Expected HTML not found: ", html_out)
if (!dir.exists("docs")) dir.create("docs", recursive = TRUE)
file.copy(html_out, file.path("docs", basename(html_out)), overwrite = TRUE)
cat("Copied to docs/\n")

cat("\n== Next: git status/commit/push checklist ==\n")
cat("1) git status\n")
cat("2) Inspect report locally:\n   open docs/fishery_sel_form_comparison.html\n")
cat("3) Stage changes (example):\n   git add -f analysis/output/fishery_sel_forms/*.rds analysis/output/fishery_sel_forms/summary.csv \\n             reporting/fishery_sel_form_comparison.qmd reporting/fishery_sel_form_comparison.html \\n             docs/fishery_sel_form_comparison.html\n")
cat("4) Commit: git commit -m \"Update scenario results + report\"\n")
cat("5) Push:   git push origin main\n")

cat("\nAll done.\n")
