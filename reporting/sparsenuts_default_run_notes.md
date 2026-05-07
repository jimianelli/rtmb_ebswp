# SparseNUTS Default Run Notes

Date: 2026-04-26

Current RTMB configuration: steepness is fixed at `0.67`.

The RTMB-ADMB SparseNUTS run was created by rendering:

```bash
quarto render reporting/ebs_pollock_rtmb_ebswp_assessment.qmd -P run_sparsenuts:true
```

The report calls the library default:

```r
SparseNUTS::sample_snuts(obj)
```

## Serial/Parallel

The call was not configured as serial. `SparseNUTS::sample_snuts()` defaults to
`chains = 4` and `cores = chains`. The saved object has four chains:

- samples dimension: `1150 x 4 x 1339`
- warmup: `150`
- post-warmup samples per chain: `1000`
- algorithm: `SNUTS`
- metric: `dense`

Each chain is still internally sequential, but the run was not requested with
`cores = 1`.

## Saved Output

- SparseNUTS RDS:
  `analysis/output/sparsenuts/rtmb_ebswp_sparsenuts_default.rds`
- Rendered report:
  `reporting/ebs_pollock_rtmb_ebswp_assessment.html`

## Diagnostics

Summary from the saved object after fixing steepness at `0.67`:

- maximum Rhat: `1.009443`
- minimum bulk ESS: `526.8279`
- minimum tail ESS: `183.0811`
- percent divergent: `0.68`
- percent treedepth: `0`
- number below diagnostic threshold: `0`

The slowest parameters by Rhat/ESS were headed by:

| Parameter | Rhat | Bulk ESS | Tail ESS |
|---|---:|---:|---:|
| `sel_a50_bts_dev[8]` | 1.009443 | 526.8279 | 183.0811 |
| `sel_slp_bts_dev[8]` | 1.008453 | 559.7373 | 187.1109 |

Treat this default run as a computational diagnostic, not accepted posterior
inference.

## Completed MCMC Plots

The report now includes the following default-run plots:

- Slow-order pairs-style plot:
  `analysis/output/figures/rtmb_ebswp_sparsenuts_pairs_slow.png`
- `SparseNUTS::plot_marginals(..., order = "slow")`:
  `analysis/output/figures/rtmb_ebswp_sparsenuts_marginals_slow.png`
- `SparseNUTS::plot_sampler_params()`:
  `analysis/output/figures/rtmb_ebswp_sparsenuts_sampler_params.png`
- `SparseNUTS::plot_Q()`:
  `analysis/output/figures/rtmb_ebswp_sparsenuts_Q.png`
- `SparseNUTS::plot_uncertainties()`:
  `analysis/output/figures/rtmb_ebswp_sparsenuts_uncertainties.png`

No installed namespace in this workspace exposes a callable `pairs_admb()`
function. The report therefore includes an equivalent slow-order pairs-style
plot from the saved SparseNUTS samples and labels that limitation explicitly.
