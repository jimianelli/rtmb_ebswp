# September 2025 ADMB-to-RTMB bridge

This branch restores the fixed-parameter EBS pollock bridge published in
September 2025. It is a numerical reference configuration, not the branch for
new selectivity, profiling, or MCMC development.

## Sources

- Published report: <https://noaa-afsc.github.io/EBS_pollock/doc/Sept_2025.html>
- NOAA-AFSC repository: <https://github.com/NOAA-AFSC/EBS_pollock>
- RTMB bridge commit: `8af13c6417429fba9725b0a0ce6bc351ba4a620d`
- Report-era source commit: `e848196`
- Final bundled ADMB template commit: `0f1a9b7a76f30d77b3970bb2c574cae7de465d53`

The restored `R/Rpm.R`, `R/model_funs.R`, `R/utilities.R`, and `R/config.R`
come from the report-era NOAA-AFSC source. The two changes to `R/config.R` are
case-only path corrections from `rpm.R` to the tracked filename `Rpm.R` so the
configuration also works on case-sensitive filesystems.

The repository root symlink `runs -> admb/runs` and the bridge alias
`admb/runs/rtmb -> for_rtmb` preserve the paths used by the historical code
without duplicating the bundled inputs.

## Artifact identity

The report-era `runs/rtmb/pm.tpl` and the bundled
`admb/runs/for_rtmb/pm.tpl` have the same SHA-256 digest:

```text
d1384a696de864c0918542dcb5992272ba3283b617f00ad8efc7b3e2a391d66b
```

The successful historical `R/Rpm.R` has SHA-256 digest:

```text
b3f360b324cc57af8248d97dea12e9f6cf94e5e07fae5494401a369e4de775cb
```

## Regression gate

Run:

```bash
Rscript tests/test_sept_2025_bridge.R
```

The gate requires:

- absolute total NLL difference no greater than `0.005`;
- maximum percent difference among key outputs no greater than `0.0006%`;
- maximum absolute RTMB gradient no greater than `1e-5`; and
- `all.equal(..., tolerance = 1e-5)` for every key output.

The published and locally reproduced benchmarks are:

```text
Absolute total NLL difference: 0.004098
Maximum absolute gradient:     6.36535e-06
Maximum key percent difference approximately 0.0005%
```

To save the complete comparison table during a run:

```bash
RTMB_BRIDGE_COMPARISON_CSV=analysis/output/sept_2025_bridge.csv \
  Rscript tests/test_sept_2025_bridge.R
```

## Development rule

Run this gate before and after each later model change. Changes to parameter
mapping, fishery selectivity, fishing mortality, state propagation, likelihood
constants, or fixed values must not be combined until the first change passes
the bridge independently.

Grant Adams' `m23_rceattle_full` is a later and intentionally different ADMB
configuration. It should be implemented as a separately named RTMB model after
this reference bridge remains closed.
