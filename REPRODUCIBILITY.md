# Reproducibility

This repository is designed to reproduce the RTMB side of the EBS pollock
RTMB-ADMB bridge using the bundled minimum ADMB bridge inputs. A different ADMB
bridge run can still be supplied through `POLLOCK_ROOT` or `POLLOCK_BASE`.

## 1. Bridge Data Layout

The repository includes this minimum layout:

```text
admb/
  runs/
    for_rtmb/
      pm.dat
      pm.par
      pm.rep
      pm.std
      pm.tpl
      control.dat
      proj_df.csv
      proj/
        tacpar.dat
    data/
      cov_2024.dat
      pm_24.dat
      selvar24.dat
      surv.dat
      surveycpue.dat
      wtage2024.dat
```

To override the bundled bridge inputs with another pollock workspace, set:

```bash
export POLLOCK_ROOT=/path/to/pollock
```

If neither `POLLOCK_ROOT` nor `POLLOCK_BASE` is set, the code uses the bundled
`admb/runs/` tree.

## 2. R Dependencies

Core model runs require:

- R
- RTMB
- tidyverse
- patchwork
- ebswp

Optional diagnostics and reports use:

- Quarto
- gt
- ggthemes
- ggridges
- cowplot
- SparseNUTS
- afscOSA
- spmR (>= 0.3.0; required for the validated seven-scenario Tier 3 table)
- generics
- tibble
- dplyr
- yardstick (optional common fit metrics)
- hardhat (optional weighted fit metrics)
- reactable
- htmlwidgets
- webshot2

Install package versions consistent with the assessment workspace used to create
the ADMB bridge files. When making archival runs, record the output of:

```r
sessionInfo()
```

or use `renv::snapshot()` in a private analysis clone to create a lockfile.

## 3. Rebuild Standard RTMB Output

From the repository root:

```bash
Rscript R/write_output.R
```

Expected output:

```text
analysis/output/base.rds
```

This file is generated and should not be committed.

## 4. Optional Diagnostics

Build standardized `tidy()`, `glance()`, and `augment()` outputs for the saved
custom RTMB, Rceattle, and SPoRC pollock fits with:

```bash
Rscript scripts/build_tidy_pollock_outputs.R
```

This cross-repository step expects the sibling directories `ebswp_rceattle`
and `sporc_ebswp` beside this repository. When `yardstick` is installed, the
script also writes grouped catch and index metrics.

Run the BTS-only comparison:

```bash
Rscript R/run_only_bts.R
```

Prepare OSA inputs and run OSA diagnostics:

```bash
Rscript R/prepare_osa_inputs.R
Rscript R/run_osa_comps.R
```

Run a five-peel retrospective:

```bash
Rscript R/run_retrospective.R
```

Run the hierarchical Form-2 SparseNUTS diagnostic with the sparse metric and
`adapt_delta = 0.95`:

```bash
SPARSENUTS_FORMS=2 SPARSENUTS_METRIC=sparse \
SPARSENUTS_ADAPT_DELTA=0.95 SPARSENUTS_OUTPUT_TAG=sparse_adapt095 \
Rscript R/run_sparsenuts_fishery_sel_forms.R
```

## 5. Render Report

After creating the required generated outputs:

```bash
quarto render
```

To force the SparseNUTS run during render:

```bash
quarto render -P run_sparsenuts:true
```

The project writes publishable HTML directly to `docs/`. Commit those two
published pages when releasing an update. HTML under `reporting/` is a local
artifact and remains ignored.

## 6. Archival Checklist

For a reproducible release, record:

- Git commit hash of this repository.
- Provenance of the bundled or external ADMB bridge workspace.
- Checksums of `pm.par`, `pm.rep`, and `pm.tpl`.
- R version and package versions.
- Commands run and parameter overrides.
- Generated output checksums for any archived RDS, HTML, or projection files.
