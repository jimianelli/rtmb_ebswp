#!/usr/bin/env Rscript

# Plot simulated fishery proportions-at-age from base-model pseudo-data.
#
# Usage:
#   NSIM=30 Rscript R/plot_simulated_fishery_age_props.R
#
# Output:
#   analysis/output/figures/simulated_fishery_age_props_by_year_30.png

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
  file.path(
    rtmb_root,
    "analysis", "output", "figures",
    sprintf("simulated_fishery_age_props_by_year_%d.png", nsim_plot)
  )
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

years <- as.numeric(base$report$yrs_fsh_data)
ages <- seq_len(ncol(base$report$phat_fsh))
expected <- as.matrix(base$report$phat_fsh)

pollock_root <- Sys.getenv("POLLOCK_ROOT", unset = NA_character_)
if (is.na(pollock_root) || !nzchar(pollock_root)) {
  pollock_root <- Sys.getenv("POLLOCK_BASE", unset = dirname(rtmb_root))
}
pollock_root <- normalizePath(pollock_root, mustWork = TRUE)
source(file.path(rtmb_root, "R", "utilities.R"))
input_data <- Get_Data()
observed <- as.matrix(input_data$oac_fsh)
if (length(years) == 0 || any(!is.finite(years))) {
  years <- as.numeric(input_data$yrs_fsh_data)
}

sim_array <- array(
  NA_real_,
  dim = c(length(years), length(ages), nsim_plot),
  dimnames = list(year = years, age = ages, sim = seq_len(nsim_plot))
)
for (i in seq_along(sim_list)) {
  sim_array[, , i] <- as.matrix(sim_list[[i]]$oac_fsh)
}

plot_years <- sort(years)
n_panel <- length(plot_years)
nrow_panel <- 10L
ncol_panel <- ceiling(n_panel / nrow_panel)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
png(out_path, width = 600 * ncol_panel, height = 360 * nrow_panel, res = 180)
op <- par(
  mfcol = c(nrow_panel, ncol_panel),
  mar = c(1.25, 1.6, 1.45, 0.35),
  oma = c(3.2, 3.8, 2.8, 0.2),
  mgp = c(1.1, 0.35, 0),
  tcl = -0.2,
  las = 1
)
on.exit(par(op), add = TRUE)

sim_col <- grDevices::adjustcolor("#4f8fc0", alpha.f = 0.24)
fit_col <- "#c43c39"
obs_col <- "#202020"

for (yr in plot_years) {
  iyear <- match(yr, years)
  y <- sim_array[iyear, , , drop = FALSE]
  y_range <- range(c(y, expected[iyear, ], observed[iyear, ]), finite = TRUE)
  plot(
    ages,
    expected[iyear, ],
    type = "n",
    ylim = y_range,
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    main = yr,
    cex.main = 0.78
  )
  matlines(ages, y[1, , ], lty = 1, lwd = 0.8, col = sim_col)
  lines(ages, expected[iyear, ], lwd = 1.8, col = fit_col)
  points(ages, observed[iyear, ], pch = 16, cex = 0.28, col = obs_col)
  axis(1, at = c(1, 5, 10, 15), labels = c(1, 5, 10, 15), cex.axis = 0.58)
  axis(2, cex.axis = 0.52)
  if (identical(yr, plot_years[1])) {
    legend(
      "topright",
      legend = c("Simulated", "Fitted", "Observed"),
      col = c(grDevices::adjustcolor("#4f8fc0", alpha.f = 0.55), fit_col, obs_col),
      lty = c(1, 1, NA),
      pch = c(NA, NA, 16),
      lwd = c(1.2, 2.1, NA),
      cex = 0.55,
      bty = "n"
    )
  }
  box()
}

mtext("Age", side = 1, outer = TRUE, line = 2.1)
mtext(
  sprintf("Base-Model Simulated Fishery Proportions at Age by Year (%d draws)", nsim_plot),
  side = 3,
  outer = TRUE,
  line = 0.7,
  cex = 1.15,
  font = 2
)
dev.off()

cat("Wrote fishery age-proportion simulation plot:", out_path, "\n")
