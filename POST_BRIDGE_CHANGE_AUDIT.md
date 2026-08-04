# Post-bridge change audit

This audit identifies when the numerical equivalence documented in the
September 2025 report changed in the standalone `rtmb_ebswp` history.

## Benchmarks

The same bundled ADMB `pm.par`, `pm.rep`, and data were evaluated at selected
source revisions. Percent differences are maximum absolute differences relative
to the ADMB report.

| Source state | NLL difference | Maximum gradient | SSB % | N % | F % | Fishery selectivity % | Catch % |
|---|---:|---:|---:|---:|---:|---:|---:|
| September 2025 archived bridge | -0.004098 | 0.00000637 | 0.000447 | 0.000486 | 0.000482 | 0.000493 | 0.000383 |
| NOAA-AFSC `215efa4` (2026-01-28) | -0.004098 | 0.0000149 | 0.000447 | 0.000486 | 0.000482 | 0.000493 | 0.000383 |
| Standalone `567e2ee` (2026-05-07) | 0.506449 | 4.68066 | 0.000447 | 0.000486 | 0.000482 | 0.000493 | 0.000383 |
| Standalone `7301aa0` (selectivity alternatives) | 0.506449 | 4.68066 | 0.000447 | 0.000486 | 0.000482 | 0.000493 | 0.000383 |
| Standalone `af5c490` | 0.506449 | 4.68066 | 0.000447 | 0.000486 | 0.000482 | 0.000493 | 0.000383 |
| Standalone `7db7217` (tie old ages) | 10.7369 | 113.409 | 3.47938 | 12.0565 | 26.5334 | 26.5336 | 6.51986 |
| Standalone `cd46353` | 10.7369 | 113.409 | 3.47938 | 12.0565 | 26.5334 | 26.5336 | 6.51986 |

## Findings

There are two separate regressions.

### Objective and gradient regression at standalone creation

The initial standalone repository preserved the ADMB state trajectories but
overrode the ADMB steepness estimate with `0.67` and changed the mapping. This
increased the objective discrepancy from approximately `0.0041` to `0.5064`
and the maximum gradient from near zero to `4.68`.

This should be handled by separating two configurations:

- `bridge`: use every ADMB parameter and the historical map unchanged;
- `assessment-development`: allow fixed steepness and revised parameter maps.

### State and prediction regression at `7db7217`

Commit `7db7217e4016c2f9d8f6b988fb308ca456b753a8`, “Tie old-age fishery
selectivity and update comparison,” is the first tested revision with large
trajectory differences. It added `cap_old_age_log_selectivity()` and applied it
to fishery selectivity form 0, the historical bridge form. Ages 11 and older
were set equal to age 11 and the curve was recentered. The corresponding ADMB
template was not changed, so F, numbers, catch, SSB, and likelihoods diverged.

The safe implementation is to leave form 0 unchanged and apply the old-age cap
only to explicitly named sensitivity forms. If a capped base model is desired,
it needs a matching ADMB reference run and must not be called the September
2025 bridge.

## Reapplication order

Starting from this restored branch:

1. Preserve the historical parameter map and pass the bridge gate.
2. Add selectable fishery forms without changing form 0; run the gate.
3. Add double-logistic and 2D-AR1 sensitivities; run the gate after each.
4. Add old-age tying only as a separate form or option; keep it off in bridge mode.
5. Add profiling helpers; run the gate.
6. Add SparseNUTS support; run the deterministic gate before sampling.
7. Add Grant's `m23_rceattle_full` as a separate configuration with its own ADMB artifacts and bridge thresholds.

The current dirty `agent/rtmbprof-stage1` branch remains untouched and should
be reconciled only after these deterministic changes have been reintroduced
behind explicit configuration switches.
