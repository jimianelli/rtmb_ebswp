#!/usr/bin/env Rscript

# Plot simulated BTS index series from the base-model pseudo-data.
#
# Usage:
#   NSIM=30 Rscript R/plot_simulated_bts_index.R
#
# Output:
#   analysis/output/figures/simulated_bts_index_30.png

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

rtmb_root <- normalizePath(getwd(), mustWork = TRUE)
nsim_plot <- as.integer(Sys.getenv("NSIM", "30"))
sim_path <- Sys.getenv(
  "SIM_FILE",
  file.path(rtmb_root, "analysis", "output", "simulated", "base_simulated_datasets.rds")
)
base_path <- Sys.getenv(
  "BASE_FILE",
  file.path(rtmb_root, "analysis", "output", "base.rds")
)
out_path <- Sys.getenv(
  "OUTPUT",
  file.path(rtmb_root, "analysis", "output", "figures", sprintf("simulated_bts_index_%d.png", nsim_plot))
)

if (!file.exists(sim_path)) {
  stop("Missing simulated datasets: ", sim_path, "\nRun: NSIM=", nsim_plot, " Rscript R/simulate_base_datasets.R")
}
if (!file.exists(base_path)) {
  stop("Missing base model output: ", base_path)
}

sims <- readRDS(sim_path)
base <- readRDS(base_path)
sim_list <- sims$simulations
if (length(sim_list) < nsim_plot) {
  stop(
    "Simulation file has ", length(sim_list), " datasets but ", nsim_plot,
    " are requested. Regenerate with: NSIM=", nsim_plot,
    " Rscript R/simulate_base_datasets.R"
  )
}
sim_list <- sim_list[seq_len(nsim_plot)]

years <- sim_list[[1]]$yrs_bts_data
sim_matrix <- vapply(sim_list, function(x) as.numeric(x$ob_bts), numeric(length(years)))
fitted_bts <- as.numeric(base$report$eb_bts)
observed_bts <- as.numeric(sim_list[[1]]$ob_bts)

# Recover the original observed BTS series from input files when available.
pollock_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
if (is.na(pollock_root) || !nzchar(pollock_root)) {
  pollock_root <- Sys.getenv("POLLOCK_BASE", unset = dirname(rtmb_root))
}
pollock_root <- normalizePath(pollock_root, mustWork = TRUE)
source(file.path(rtmb_root, "R", "utilities.R"))
input_data <- try(Get_Data(), silent = TRUE)
if (!inherits(input_data, "try-error") && !is.null(input_data$ob_bts)) {
  observed_bts <- as.numeric(input_data$ob_bts)
}

ylim <- range(c(sim_matrix, fitted_bts, observed_bts), finite = TRUE)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
png(out_path, width = 1800, height = 1100, res = 180)
par(mar = c(4.5, 4.8, 2.2, 1), las = 1)
plot(
  years,
  fitted_bts,
  type = "n",
  ylim = ylim,
  xlab = "Year",
  ylab = "BTS biomass index",
  main = sprintf("Base-Model Simulated BTS Index (%d draws)", nsim_plot)
)
matlines(years, sim_matrix, lty = 1, lwd = 1.1, col = grDevices::adjustcolor("#4f8fc0", alpha.f = 0.28))
lines(years, fitted_bts, lwd = 2.4, col = "#c43c39")
points(years, observed_bts, pch = 16, cex = 0.72, col = "#202020")
legend(
  "topleft",
  legend = c("Simulated datasets", "Fitted base expectation", "Observed input"),
  col = c(grDevices::adjustcolor("#4f8fc0", alpha.f = 0.5), "#c43c39", "#202020"),
  lty = c(1, 1, NA),
  pch = c(NA, NA, 16),
  lwd = c(1.2, 2.4, NA),
  bty = "n"
)
box()
dev.off()

cat("Wrote BTS simulation plot:", out_path, "\n")
