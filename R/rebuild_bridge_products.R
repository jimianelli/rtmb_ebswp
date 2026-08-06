#!/usr/bin/env Rscript

# Rebuild all published products from the September 2025 ADMB-to-RTMB bridge.
# Run from the repository root:
#   Rscript R/rebuild_bridge_products.R

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "R", "config.R"))) {
  stop("Run this script from the rtmb_ebswp repository root.")
}

run_r <- function(script) {
  message("\n==> Rscript ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) stop(script, " failed with status ", status)
}

run_r(file.path("R", "write_output.R"))
run_r(file.path("tests", "test_sept_2025_bridge.R"))
run_r(file.path("R", "run_retrospective.R"))
run_r(file.path("R", "run_spmr_tier3_projection.R"))
run_r(file.path("R", "run_sparsenuts_default.R"))

base_file <- file.path("analysis", "output", "base.rds")
base_md5 <- unname(tools::md5sum(base_file))
retro_md5 <- readRDS(file.path("analysis", "output", "retro_9_peel.rds"))$
  base_lineage$md5
projection_md5 <- read.csv(
  file.path("analysis", "output", "spmR_projection", "base_lineage.csv")
)$model_md5[1]
alt2_md5 <- read.csv(
  file.path(
    "analysis", "output", "spmR_projection_alt2_fixed1300",
    "base_lineage.csv"
  )
)$model_md5[1]
mcmc_md5 <- read.csv(
  file.path("reporting", "data", "sparsenuts_base_run_summary.csv")
)$base_md5[1]

lineage <- c(
  base = base_md5,
  retrospective = retro_md5,
  projections = projection_md5,
  alternative_2 = alt2_md5,
  SparseNUTS = mcmc_md5
)
if (!all(lineage == base_md5)) {
  stop("Downstream base-file lineage mismatch:\n", paste(names(lineage), lineage, collapse = "\n"))
}

message("\n==> quarto render reporting/ebs_pollock_rtmb_ebswp_assessment.qmd")
quarto_status <- system2(
  "quarto",
  c("render", file.path("reporting", "ebs_pollock_rtmb_ebswp_assessment.qmd"))
)
if (!identical(quarto_status, 0L)) stop("Quarto render failed.")

rendered <- file.path("reporting", "ebs_pollock_rtmb_ebswp_assessment.html")
published <- file.path("docs", "ebs_pollock_rtmb_ebswp_assessment.html")
if (!file.copy(rendered, published, overwrite = TRUE)) {
  stop("Could not copy rendered assessment to docs/.")
}

message("\nRebuild complete. Common base checksum: ", base_md5)

