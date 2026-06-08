# Data Availability

This repository includes the minimal ADMB bridge bundle needed for standard RTMB
bridge runs. The default in-repository layout is:

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

To run against a different ADMB bridge, provide an external pollock workspace and
set `POLLOCK_ROOT` to that directory. `POLLOCK_BASE` is also recognized for
older scripts.

Generated outputs from this repository are written under `analysis/output/` and
are ignored by git.
