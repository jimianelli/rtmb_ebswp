# SparseNUTS base-model run notes

The current MCMC result was generated from the accepted base RTMB model
(`fishery_sel_form = 0`) with:

```r
SparseNUTS::sample_snuts(
  obj,
  globals = list(data = model_data),
  model_name = "EBS pollock base RTMB",
  control = list(adapt_delta = 0.95)
)
```

The target acceptance probability is increased to 0.95 to reduce divergent
transitions; other sampler settings remain at package defaults. The `globals` argument exports the RTMB data
object required to rebuild the objective in the package-default parallel worker
sessions. Run the analysis from the repository root with:

```sh
Rscript R/run_sparsenuts_default.R
quarto render reporting/ebs_pollock_rtmb_ebswp_assessment.qmd
cp reporting/ebs_pollock_rtmb_ebswp_assessment.html docs/
```

## Current run

- SparseNUTS 1.0.2, GitHub commit
  `2f3f1626219afce68fa2da0d884d4f2dca138117`
- automatic dense metric
- four chains running on four parallel workers
- 150 warmup and 1,000 retained iterations per chain
- `adapt_delta = 0.95`
- maximum R-hat: 1.006
- minimum bulk effective sample size: 1,569
- 7 warmup divergences and 1 divergence among 4,000 retained transitions

The increased target acceptance probability reduced the post-warmup
divergences from 30 in the original untuned run to one in this rerun. The
remaining transition should be located and substantive posterior checks and
comparison with deterministic uncertainty remain necessary before final use.

## Outputs

- full local fit:
  `analysis/output/sparsenuts/rtmb_ebswp_sparsenuts_default.rds`
- tracked run summary:
  `reporting/data/sparsenuts_base_run_summary.csv`
- tracked parameter monitor:
  `reporting/data/sparsenuts_base_monitor.csv`

The Quarto report reads the full fit and regenerates every MCMC table and figure
on each render, so a new completed run cannot leave old MCMC graphics in the
published assessment page.
