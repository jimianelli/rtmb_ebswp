args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
rtmb_root <- if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  normalizePath(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
osa_input_file <- Sys.getenv(
  "RTMB_OSA_INPUT_FILE",
  unset = file.path(rtmb_root, "analysis", "output", "osa_inputs.rds")
)
osa_output_dir <- Sys.getenv(
  "RTMB_OSA_OUTPUT_DIR",
  unset = file.path(rtmb_root, "analysis", "output", "osa")
)
osa_output_file <- Sys.getenv(
  "RTMB_OSA_OUTPUT_FILE",
  unset = file.path(osa_output_dir, "rtmb_ebswp_osa_residuals.rds")
)

if (!file.exists(osa_input_file)) {
  stop("Missing OSA input file at: ", osa_input_file)
}

suppressPackageStartupMessages({
  library(afscOSA)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

osa_inputs <- readRDS(osa_input_file)
input_lineage <- attr(osa_inputs, "lineage")
dir.create(osa_output_dir, recursive = TRUE, showWarnings = FALSE)

osa_runs <- lapply(names(osa_inputs), function(fleet_name) {
  input <- osa_inputs[[fleet_name]]
  afscOSA::run_osa(
    obs = input$obs,
    exp = input$exp,
    N = input$n_eff,
    fleet = fleet_name,
    index = input$index,
    years = input$years,
    index_label = "Age"
  )
})
names(osa_runs) <- names(osa_inputs)

osa_plot_components <- afscOSA::plot_osa(
  osa_runs,
  plot = FALSE,
  outpath = NULL,
  figheight = 9,
  figwidth = 12
)

osa_plot_components <- lapply(
  osa_plot_components,
  function(plot) plot + ggthemes::theme_few()
)

osa_plots <- cowplot::plot_grid(
  osa_plot_components$aggcomp,
  osa_plot_components$qq,
  osa_plot_components$bubble + ggplot2::theme(legend.position = "top"),
  osa_plot_components$bubble_pearson + ggplot2::theme(legend.position = "none"),
  nrow = 4,
  rel_heights = c(5, 5, 4.5, 4)
)

ggplot2::ggsave(
  filename = file.path(osa_output_dir, "osa_age_diagnostics.png"),
  plot = osa_plots,
  units = "in",
  bg = "white",
  height = 9,
  width = 12,
  dpi = 300
)

osa_summary <- dplyr::bind_rows(lapply(osa_runs, function(x) x$res)) |>
  dplyr::group_by(fleet) |>
  dplyr::summarise(
    residuals = dplyr::n(),
    mean = mean(resid, na.rm = TRUE),
    sdnr = stats::sd(resid, na.rm = TRUE),
    lower_2.5_pct = stats::quantile(resid, 0.025, na.rm = TRUE),
    upper_97.5_pct = stats::quantile(resid, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
utils::write.csv(
  osa_summary,
  file.path(osa_output_dir, "osa_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    runs = osa_runs,
    plots = osa_plots,
    summary = osa_summary,
    created = Sys.time(),
    r_version = R.version.string,
    afscOSA_version = as.character(utils::packageVersion("afscOSA")),
    input_md5 = unname(tools::md5sum(osa_input_file)),
    lineage = input_lineage
  ),
  osa_output_file
)

message("Wrote OSA residuals to ", osa_output_file)
