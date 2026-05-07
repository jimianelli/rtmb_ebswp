# Reproducibility

This repository is designed to reproduce the RTMB side of the EBS pollock
RTMB-ADMB bridge when supplied with the external ADMB bridge inputs.

## 1. External Data Layout

Create or identify a pollock workspace outside this repository with this
minimum layout:

```text
pollock/
  admb/
    runs/
      for_rtmb/
        pm.par
        pm.rep
        pm.tpl
      data/
```

Then set:

```bash
export POLLOCK_ROOT=/path/to/pollock
```

The repository does not version these files because they are assessment inputs
and generated ADMB bridge outputs.

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
- spmR
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

Run the default SparseNUTS diagnostic:

```bash
Rscript R/run_sparsenuts_default.R
```

## 5. Render Report

After creating the required generated outputs:

```bash
quarto render reporting/ebs_pollock_rtmb_admb_assessment.qmd
```

To force the SparseNUTS run during render:

```bash
quarto render reporting/ebs_pollock_rtmb_admb_assessment.qmd -P run_sparsenuts:true
```

Rendered HTML is generated output and should not be committed.

## 6. Archival Checklist

For a reproducible release, record:

- Git commit hash of this repository.
- Location and provenance of the external pollock workspace.
- Checksums of `pm.par`, `pm.rep`, and `pm.tpl`.
- R version and package versions.
- Commands run and parameter overrides.
- Generated output checksums for any archived RDS, HTML, or projection files.
