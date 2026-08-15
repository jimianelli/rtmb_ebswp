#!/bin/sh
set -eu

for report_name in \
  ebs_pollock_rtmb_ebswp_assessment \
  appendix_fishery_selectivity
do
  staged_file="docs/reporting/${report_name}.html"
  published_file="docs/${report_name}.html"
  if [ -f "$staged_file" ]; then
    mv "$staged_file" "$published_file"
  fi
done

rmdir docs/reporting 2>/dev/null || true

# Publish stable, reader-facing CSV products referenced by the assessment.
mkdir -p docs/data-output
cp analysis/output/tidy/model_glance.csv \
  docs/data-output/model_glance.csv
cp analysis/output/tidy/model_parameters.csv \
  docs/data-output/model_parameters.csv
cp analysis/output/tidy/model_observations.csv \
  docs/data-output/model_observations.csv
cp analysis/output/bts_age_data_bridge/comparison.csv \
  docs/data-output/bts_age_data_bridge_comparison.csv
cp analysis/output/bts_age_data_bridge/diagnostics.csv \
  docs/data-output/bts_age_data_bridge_diagnostics.csv
cp analysis/output/bts_age_data_bridge/timeseries.csv \
  docs/data-output/bts_age_data_bridge_timeseries.csv
cp analysis/output/bts_age_data_bridge/tinyvast_retro_diagnostics.csv \
  docs/data-output/tinyvast_retro_diagnostics.csv
cp analysis/output/bts_age_data_bridge/tinyvast_retro_mohn.csv \
  docs/data-output/tinyvast_retro_mohn.csv
cp analysis/output/corrected_full_age_bts/downstream_lineage.csv \
  docs/data-output/downstream_lineage.csv
cp analysis/output/corrected_full_age_bts/fishery_sel_forms/summary.csv \
  docs/data-output/fishery_selectivity_summary.csv
cp analysis/output/corrected_full_age_bts/osa/osa_summary.csv \
  docs/data-output/osa_summary.csv
cp analysis/output/corrected_full_age_bts/retrospective_data_availability.csv \
  docs/data-output/retrospective_data_availability.csv
cp analysis/output/corrected_full_age_bts/spmR_projection/spm_summary.csv \
  docs/data-output/spmr_projection_summary.csv
cp analysis/output/corrected_full_age_bts/spmR_projection/tier3_seven_scenario_table.csv \
  docs/data-output/tier3_seven_scenario_table.csv
cp reporting/data/lf_length_frequency_summary.csv \
  docs/data-output/lf_length_frequency_summary.csv
