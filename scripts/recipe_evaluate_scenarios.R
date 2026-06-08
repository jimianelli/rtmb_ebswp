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
#   Rscript scripts/recipe_evaluate_scenarios.R
#   POLLOCK_ROOT=/path/to/pollock Rscript scripts/recipe_evaluate_scenarios.R
#
# Notes:
#   - This assumes you have already set up rtmb_ebswp dependencies.
#   - This script does not modify model code; it only runs and validates outputs.

preferred_rscript <- Sys.getenv(
  "RTMB_RSCRIPT",
  unset = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
)
if (file.exists(preferred_rscript) &&
    normalizePath(R.home(), mustWork = TRUE) !=
      normalizePath("/Library/Frameworks/R.framework/Resources", mustWork = TRUE) &&
    Sys.getenv("RTMB_REEXECED", unset = "0") != "1") {
  this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  if (is.na(this_file) || !nzchar(this_file)) {
    stop("Cannot re-exec recipe because the script path was not available.")
  }
  env <- c(
    paste0("RTMB_REEXECED=1"),
    paste0("RTMB_RSCRIPT=", preferred_rscript),
    paste0("POLLOCK_ROOT=", Sys.getenv("POLLOCK_ROOT", unset = ""))
  )
  status <- system2(preferred_rscript, c(normalizePath(this_file, mustWork = TRUE), commandArgs(TRUE)), env = env)
  quit(save = "no", status = status)
}

rtmb_rscript <- if (file.exists(preferred_rscript)) preferred_rscript else file.path(R.home("bin"), "Rscript")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# helper
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

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

resolve_pollock_root <- function(repo_root) {
  has_bridge <- function(path) {
    if (is.na(path) || !nzchar(path)) return(FALSE)
    file.exists(file.path(path, "admb", "runs", "for_rtmb", "pm.par")) &&
      file.exists(file.path(path, "admb", "runs", "for_rtmb", "pm.rep")) &&
      file.exists(file.path(path, "admb", "runs", "for_rtmb", "pm.tpl")) &&
      file.exists(file.path(path, "admb", "runs", "data", "pm_24.dat"))
  }

  env_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
  if (is.na(env_root) || !nzchar(env_root)) {
    env_root <- Sys.getenv("POLLOCK_BASE", unset = NA_character_)
  }
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  candidates <- c(repo_root, dirname(repo_root), file.path(dirname(repo_root), "pollock"))
  for (cand in candidates) {
    if (has_bridge(cand)) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }

  stop("Cannot locate pollock bridge inputs. Set POLLOCK_ROOT or use bundled admb/runs files.")
}

pollock_root <- resolve_pollock_root(repo_root)
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
cmd <- sprintf(
  "POLLOCK_ROOT=%s %s R/run_fishery_selectivity_forms.R",
  shQuote(pollock_root),
  shQuote(rtmb_rscript)
)
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
    form = as.integer(forms[i]),
    path = paths[i],
    label = x$label %||% paste0("form ", forms[i]),
    objective = x$objective,
    convergence = x$convergence,
    max_gradient = x$max_gradient,
    seconds = x$seconds,
    report = x$report
  )
})

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
  if (is.null(mat)) return(tibble(component = name, ok = NA, min = NA, max = NA, max_abs_dev = NA))
  rs <- rowSums(mat)
  tibble(
    component = name,
    ok = all(is.finite(rs)) && max(abs(rs - 1)) < 1e-6,
    min = min(rs),
    max = max(rs),
    max_abs_dev = max(abs(rs - 1))
  )
}

agecomp_diag <- bind_rows(lapply(runs, function(r) {
  rep <- r$report
  bind_rows(
    # Observed
    check_comp(rep$oac_fsh, "oac_fsh"),
    check_comp(rep$oac_bts, "oac_bts"),
    check_comp(rep$oac_ats, "oac_ats"),
    # Expected
    check_comp(rep$phat_fsh, "phat_fsh"),
    check_comp(rep$phat_bts, "phat_bts"),
    check_comp(rep$phat_ats, "phat_ats")
  ) |>
    mutate(fishery_sel_form = r$form, label = r$label)
})) |>
  select(fishery_sel_form, label, component, ok, min, max, max_abs_dev)

# Print to console
for (r in runs) {
  cat("\nScenario form=", r$form, " (", r$label, ")\n", sep = "")
  print(agecomp_diag |> filter(fishery_sel_form == r$form))
}

# ---- 4) Catch fit quick check (+ save diagnostics table) ----
cat("\n== Catch fit quick check ==\n")

catch_diag <- bind_rows(lapply(runs, function(r) {
  rep <- r$report
  obs <- rep$obs_catch
  pred <- rep$pred_catch
  if (is.null(obs) || is.null(pred)) {
    return(tibble(fishery_sel_form = r$form, label = r$label,
                  corr_obs_pred = NA_real_, mean_abs_log_ratio = NA_real_))
  }
  yrs <- as.integer(names(obs) %||% names(pred) %||% seq_along(obs))
  df <- tibble(year = yrs, obs = as.numeric(obs), pred = as.numeric(pred)) |>
    filter(is.finite(obs), is.finite(pred))
  tibble(
    fishery_sel_form = r$form,
    label = r$label,
    corr_obs_pred = cor(df$obs, df$pred),
    mean_abs_log_ratio = mean(abs(log((df$obs + 1e-12) / (df$pred + 1e-12))))
  )
}))

for (i in seq_len(nrow(catch_diag))) {
  cat("form=", catch_diag$fishery_sel_form[i], ": corr(obs,pred)=", round(catch_diag$corr_obs_pred[i], 3),
      ", mean(|log(obs/pred)|)=", round(catch_diag$mean_abs_log_ratio[i], 3), "\n", sep = "")
}

# ---- 4b) Write diagnostics CSVs ----
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
write_csv(agecomp_diag, file.path(results_dir, "diagnostics_agecomp_sumto1.csv"))
write_csv(catch_diag, file.path(results_dir, "diagnostics_catch_fit.csv"))
cat("Wrote diagnostics CSVs to ", results_dir, ":\n", sep = "")
cat("- diagnostics_agecomp_sumto1.csv\n")
cat("- diagnostics_catch_fit.csv\n")

# ---- 5) Render report ----
cat("\n== Render comparison report ==\n")
qmd <- file.path("reporting", "fishery_sel_form_comparison.qmd")
if (!file.exists(qmd)) stop("Missing: ", qmd)

rtmb_r <- file.path(dirname(rtmb_rscript), "R")
quarto_home <- file.path(tempdir(), "rtmb_ebswp_quarto_home")
dir.create(quarto_home, recursive = TRUE, showWarnings = FALSE)
status <- system(sprintf(
  "HOME=%s R_LIBS_USER=%s QUARTO_R=%s quarto render %s --no-cache",
  shQuote(quarto_home),
  shQuote(Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "Library", "R", "4.6", "library"))),
  shQuote(rtmb_r),
  shQuote(qmd)
))
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
