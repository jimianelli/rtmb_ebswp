#!/bin/sh
set -eu

for report_name in \
  ebs_pollock_rtmb_ebswp_assessment \
  appendix_fishery_selectivity
do
  for extension in html pdf
  do
    staged_file="docs/reporting/${report_name}.${extension}"
    published_file="docs/${report_name}.${extension}"
    if [ -f "$staged_file" ]; then
      mv "$staged_file" "$published_file"
    fi
  done
done

if [ -f reporting/model_bridge_motherhood.html ]; then
  cp reporting/model_bridge_motherhood.html docs/index.html
  cp reporting/model_bridge_motherhood.html docs/model_bridge_motherhood.html
elif [ -f docs/reporting/model_bridge_motherhood.html ]; then
  cp docs/reporting/model_bridge_motherhood.html docs/index.html
  cp docs/reporting/model_bridge_motherhood.html docs/model_bridge_motherhood.html
elif [ -f docs/model_bridge_motherhood.html ]; then
  cp docs/model_bridge_motherhood.html docs/index.html
fi

if [ -f reporting/model_bridge_motherhood.pdf ]; then
  cp reporting/model_bridge_motherhood.pdf docs/plan_team_model_development_overview.pdf
elif [ -f docs/reporting/model_bridge_motherhood.pdf ]; then
  cp docs/reporting/model_bridge_motherhood.pdf docs/plan_team_model_development_overview.pdf
elif [ -f docs/model_bridge_motherhood.pdf ]; then
  cp docs/model_bridge_motherhood.pdf docs/plan_team_model_development_overview.pdf
fi

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
cp reporting/data/fishery_selectivity_summary.csv \
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

# Preserve the historical parent-site URL as a redirect. The report and PDF
# themselves are published only from noaa-afsc/ebswp_rceattle.
cp reporting/rceattle_redirect.html.template docs/rceattle_ebswp.html
rm -f docs/rceattle_ebswp.pdf
