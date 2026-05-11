# SMOG-CHECK Validation Notes

This directory records SMOG3-facing notes for the original SMOG-CHECK harness.
The harness resources currently remain in `SMOG-CHECK/` so the verified
drop-in validation command keeps working without path rewrites.

Run:

```bash
bash scripts/run_smogcheck_dropin_smog3.sh 1 115
```

The command creates ignored run artifacts under `smogcheck_dropin_runs/`.
