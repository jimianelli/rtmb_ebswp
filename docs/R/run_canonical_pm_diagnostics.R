#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), ".r-lib-rceattle-5.8.1"), .libPaths()))

suppressPackageStartupMessages(library(Rceattle))
stopifnot(packageVersion("Rceattle") == package_version("5.8.1"))

output_dir <- file.path("results", "canonical_pm")
method_file <- file.path(output_dir, "ebs_pollock_method_fits.rds")
fits <- readRDS(method_file)
fit <- fits$nonparametric_pm

options(mc.cores = max(1L, min(8L, parallel::detectCores() - 2L)))
set.seed(20240814)

message("Computing internal OSA residuals for the Rceattle 5.8.1 fit.")
osa <- osa_residuals(
  fit,
  source = c("index", "catch", "comp"),
  parallel = TRUE,
  seed = 20240814
)
diagnostics <- osa_diagnostics(osa)

saveRDS(osa, file.path(output_dir, "ebs_pollock_osa_residuals.rds"))
saveRDS(
  diagnostics,
  file.path(output_dir, "ebs_pollock_osa_diagnostics.rds")
)
write.csv(
  as.data.frame(diagnostics),
  file.path(output_dir, "ebs_pollock_osa_diagnostics.csv"),
  row.names = FALSE
)

message("Rceattle 5.8.1 OSA diagnostics saved under ", output_dir, ".")
