# EBS Pollock RTMB Bridge

This repository contains a standalone RTMB reimplementation of the ADMB EBS
pollock bridge model. It is intended for model-port debugging, reproducibility
checks, and diagnostics comparing RTMB output to a dedicated ADMB bridge run.

The repository contains source code, reporting scripts, and the minimal ADMB
bridge bundle needed for standard RTMB runs. Generated RDS files, rendered HTML
reports, full ADMB run products, and projection executables are intentionally
not versioned here.

## Bridge Inputs

The default run uses the bundled bridge files in this repository:

- `admb/runs/for_rtmb/pm.par`
- `admb/runs/for_rtmb/pm.rep`
- `admb/runs/for_rtmb/pm.tpl`
- `admb/runs/for_rtmb/pm.dat`
- `admb/runs/data/`

To compare against a different ADMB bridge run, point the repository at another
pollock workspace with:

```bash
export POLLOCK_ROOT=/path/to/pollock
```

For compatibility with older scripts, `POLLOCK_BASE` is also recognized when
`POLLOCK_ROOT` is not set. If neither variable is set, the code uses the
in-repository `admb/runs/` bundle.

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

## Likelihood profiles

The profile workflow fixes a selected parameter at each grid point and either
reoptimizes all remaining parameters (a likelihood profile) or evaluates the
objective without reoptimization (a slice). It records the total objective,
optimizer diagnostics, and the contribution from each reported likelihood
component.

Run the default 17-point reoptimized profile for `log_avgrec` with:

```bash
Rscript R/run_likelihood_profiles.R
```

Configure a run with environment variables. Repeated RTMB parameter names must
use one-based occurrence notation such as `log_rec_devs[10]`.

```bash
PROFILE_PARAMETERS="log_Rzero,log_q_ats" \
PROFILE_POINTS=11 \
PROFILE_HALF_WIDTH=0.3 \
PROFILE_MODE=reopt \
Rscript R/run_likelihood_profiles.R
```

When a requested parameter is fixed in `R/config.R`, the profile runner rebuilds
the RTMB map with that parameter released before fitting the profile base. This
is the default behavior for `log_avgrec`; it changes the diagnostic
configuration but does not alter the standard bridge configuration used by
other scripts.

The `log_avgrec` default spans plus or minus 2 units on the log scale. This wider
range is intentional: the recruitment deviations compensate for modest changes
in average recruitment, and the narrower plus-or-minus-0.35 trial reached only
about 0.08 Delta NLL. The wider profile crosses 1.92 on both sides.

`PROFILE_HALF_WIDTH` is measured on the fitted parameter scale, which is the log
scale for parameters whose names begin with `log_`. Outputs are written under
the ignored `analysis/output/profiles/` directory as RDS, CSV, and PNG files.
The CSV records convergence codes and maximum absolute gradients; inspect these
before interpreting profile shape. `PROFILE_MODE=slice` is useful for a quick
code check but is not a replacement for reoptimization.

Profile figures use `ggthemes::theme_few()`, include the total objective, and
use a common 0--2.1 Delta NLL scale for every facet. The Objective panel marks
Delta NLL = 1.92, the usual approximate 95% likelihood interval threshold for
one profiled parameter. Choose a grid wide enough to cross that threshold on
both sides of the minimum; increase `PROFILE_HALF_WIDTH` when it does not.

The runner requires a base-fit convergence code of zero and a maximum absolute
gradient no larger than `0.002`. Increase `PROFILE_MAX_EVAL` (default `5000`) if
the optimizer stops early. `PROFILE_GRADIENT_TOL` changes that threshold.
`PROFILE_ALLOW_NONCONVERGED=true` bypasses the check for exploratory debugging,
but output from such a run should not be used as assessment evidence.

Run the focused helper tests with:

```bash
Rscript tests/test_profile_components.R
```

### Strict RTMBprof stock-scale profile

The Stage 1 RTMBprof workflow profiles `log_q_ats`, the catchability for the
acoustic-trawl survey, as the default stock-scale diagnostic:

```bash
Rscript R/run_rtmbprof_stage1.R
```

The runner refits the default mapped model, requires a converged base fit,
checks that every named model component sums to the joint objective, profiles
the focal parameter with nuisance reoptimization, and records changes in
terminal and mean spawning biomass, unfished spawning biomass, and the
analytically derived BTS catchability. Outputs are written under the ignored
`analysis/output/rtmbprof_stage1/` directory.

The model's `log_q_bts` parameter is not a valid profile target in the current
configuration. It is fixed in the map and the effective BTS catchability is
subsequently calculated as `mean(ob_bts) / mean(eb_bts)`. Profiling it would
therefore require changing the fitted model. Short candidate profiles found
that `log_q_ats` has broader leverage on population scale than `log_q_cpue`.
`log_q_avo` has somewhat greater leverage on terminal spawning biomass but
less on mean spawning biomass and represents the shorter acoustic
vessels-of-opportunity series.

Run the strict component-closure smoke test with:

```bash
Rscript tests/test_rtmbprof_stage1.R
```

### Posterior-objective comparison

Stage 2 evaluates the same complete joint objective at the saved SparseNUTS
draws and compares its binned conditional summaries with the Stage 1
constrained `log_q_ats` profile:

```bash
Rscript R/run_rtmbprof_stage2.R
```

The runner explicitly removes the warmup iterations recorded in the saved
`tmbfit`, maps sampler columns into the fitted RTMB parameter order, and checks
component closure at every retained draw. It writes diagnostics, conditional
summaries, correlations with other scale-related parameters, and profile/MCMC
figures under `analysis/output/rtmbprof_stage2/`. The MCMC summaries retain
nuisance-parameter variation and are posterior-objective diagnostics; they are
not likelihood profiles. In this 1,338-dimensional fit, their raw objective
values are expected to lie roughly \(1338/2 = 669\) units above the mode under
a local Gaussian approximation, so the runner also produces a separately
centred shape comparison. For a quick smoke run, set
`RTMBPROF_MAX_DRAWS=40`.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the full workflow and expected
software dependencies.

The rendered report is published with GitHub Pages at:

```text
https://jimianelli.github.io/rtmb_ebswp/
```

## Repository Contents

- `R/`: RTMB model, data utilities, run scripts, and output writers.
- `analysis/`: top-level analysis entry points.
- `admb/runs/`: minimal ADMB bridge inputs and comparison outputs.
- `reporting/`: Quarto diagnostic report source and run notes.
- `data/README.md`: description of bundled bridge data and override behavior.

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
