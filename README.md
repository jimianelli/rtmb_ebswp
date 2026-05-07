# EBS Pollock RTMB Bridge

This repository contains a standalone RTMB reimplementation of the ADMB EBS
pollock bridge model. It is intended for model-port debugging, reproducibility
checks, and diagnostics comparing RTMB output to a dedicated ADMB bridge run.

The repository contains source code and reporting scripts only. Assessment
inputs, ADMB outputs, generated RDS files, rendered HTML reports, and projection
executables are intentionally not versioned here.

## External Inputs

The RTMB model is initialized and compared against an external pollock workspace
that must contain:

- `admb/runs/for_rtmb/pm.par`
- `admb/runs/for_rtmb/pm.rep`
- `admb/runs/for_rtmb/pm.tpl`
- `admb/runs/data/`

Point the repository at that workspace with:

```bash
export POLLOCK_ROOT=/path/to/pollock
```

For compatibility with older scripts, `POLLOCK_BASE` is also recognized when
`POLLOCK_ROOT` is not set.

## Run

From this repository root:

```bash
Rscript R/write_output.R
```

The standard output is written to:

```text
analysis/output/base.rds
```

To fit the model directly:

```bash
Rscript analysis/Run_rpm.R
```

To render the diagnostics report after generating required outputs:

```bash
quarto render reporting/ebs_pollock_rtmb_ebswp_assessment.qmd
```

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the full workflow and expected
software dependencies.

The rendered report is published with GitHub Pages at:

```text
https://jimianelli.github.io/rtmb_ebswp/
```

## Repository Contents

- `R/`: RTMB model, data utilities, run scripts, and output writers.
- `analysis/`: top-level analysis entry points.
- `reporting/`: Quarto diagnostic report source and run notes.
- `data/README.md`: description of external input requirements.

Generated outputs are ignored by git and should remain under
`analysis/output/` or `output/`.

## Status

This is a bridge implementation, not a final stock assessment product. Results
should not be interpreted as accepted assessment advice without the usual
convergence, estimability, retrospective, sensitivity, and review diagnostics.

## License

Code in this repository is available under the MIT License. See
[LICENSE.md](LICENSE.md). External assessment data, ADMB bridge outputs, and
generated products are not included in this license unless explicitly stated by
their source.
