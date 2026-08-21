# SparseNUTS hierarchical Form-2 run notes

The current MCMC result was generated from the hierarchical Form-2
double-logistic RTMB model (`fishery_sel_form = 2`) with a sparse metric and
`adapt_delta = 0.95`. The run used the accepted RTMB bridge output as its base
lineage; it is a selectivity sensitivity rather than the accepted base model.

```r
SparseNUTS::sample_snuts(
  obj,
  metric = "sparse",
  globals = list(data = data),
  control = list(adapt_delta = 0.95)
)
```

The target acceptance probability is increased to 0.95 to reduce divergent
transitions. The `globals` argument exports the RTMB data object required to
rebuild the objective in the parallel worker sessions. Run the analysis from
the repository root with:

```sh
SPARSENUTS_FORMS=2 SPARSENUTS_METRIC=sparse \
SPARSENUTS_ADAPT_DELTA=0.95 SPARSENUTS_OUTPUT_TAG=sparse_adapt095 \
Rscript R/run_sparsenuts_fishery_sel_forms.R
quarto render
```

## Current run

- SparseNUTS 1.0.2, GitHub commit
  `2f3f1626219afce68fa2da0d884d4f2dca138117`
- sparse metric
- four chains running on four parallel workers
- 150 warmup and 1,000 retained iterations per chain
- `adapt_delta = 0.95`
- maximum R-hat: 1.0065
- minimum bulk effective sample size: 1,199
- minimum tail effective sample size: 363
- 3 divergences among 4,000 retained transitions (0.07%)

The remaining divergent transitions should be located, and substantive
posterior checks and comparison with deterministic uncertainty remain
necessary before final use.

## Outputs

- full local fit:
  `analysis/output/sparsenuts/fishery_sel_forms/rtmb_ebswp_sparsenuts_form_2_sparse_adapt095.rds`
- tracked run summary:
  `reporting/data/sparsenuts_base_run_summary.csv`
- tracked parameter monitor:
  `reporting/data/sparsenuts_base_monitor.csv`

The Quarto report reads the full fit and regenerates every MCMC table and figure
on each render, so a new completed run cannot leave old MCMC graphics in the
published assessment page.
