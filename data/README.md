# Data Availability

This repository does not include assessment input data or ADMB bridge outputs.

To run the code, provide an external pollock workspace and set `POLLOCK_ROOT`.
The expected minimum layout is:

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

Generated outputs from this repository are written under `analysis/output/` and
are ignored by git.
