#!/usr/bin/env Rscript

source(file.path("R", "tidy-pollock-fit.R"))

pollock_root <- normalizePath("..", mustWork = TRUE)
paths <- list(
  rtmb = file.path("analysis", "output", "fishery_sel_forms", "fishery_sel_form_0.rds"),
  rceattle = file.path(pollock_root, "ebswp_rceattle", "results", "ebs_pollock_method_fits.rds"),
  sporc = file.path(pollock_root, "sporc_ebswp", "analysis", "outputs", "family_sporc.rds")
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing model artifacts: ", paste(missing, collapse = ", "))

fits <- list(
  rtmb = as_pollock_fit_rtmb(readRDS(paths$rtmb), model = "base coefficients"),
  rceattle = as_pollock_fit_rceattle(readRDS(paths$rceattle), model = "nonparametric_pm"),
  sporc = as_pollock_fit_sporc(readRDS(paths$sporc), model = "family_sporc")
)

out_dir <- file.path("analysis", "output", "tidy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(do.call(rbind, lapply(fits, function(x) as.data.frame(generics::glance(x)))),
          file.path(out_dir, "model_glance.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(fits, function(x) as.data.frame(generics::tidy(
  x, parameters = c("fixed_effect", "random_effect", "derived_quantity", "fixed_input")
)))), file.path(out_dir, "model_parameters.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(fits, function(x) as.data.frame(generics::augment(x)))),
          file.path(out_dir, "model_observations.csv"), row.names = FALSE)
saveRDS(fits, file.path(out_dir, "pollock_fits.rds"))

if (requireNamespace("yardstick", quietly = TRUE)) {
  metrics <- do.call(rbind, lapply(fits, pollock_fit_metrics))
  write.csv(as.data.frame(metrics), file.path(out_dir, "yardstick_metrics.csv"), row.names = FALSE)
} else {
  message("yardstick is unavailable; tidy outputs were written without yardstick_metrics.csv")
}

message("Wrote common pollock outputs to ", normalizePath(out_dir))
