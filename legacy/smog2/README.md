# Legacy SMOG2 Layout

The original SMOG2 source, templates, SMOG-CHECK harness, and historical
examples are intentionally preserved in their existing repository locations for
this cleanup pass.

Why not move them yet?

- The drop-in SMOG-CHECK validation command has proven paths into `SMOG-CHECK/`,
  `share/`, `src/`, and `src/tools/`.
- The two-stage parity scripts use the current layout to compare official
  Docker SMOG2 baselines against local SMOG3 outputs.
- Moving those assets safely requires a coordinated path migration across
  scripts, tests, documentation, and packaging.

SMOG3 runtime files that are needed by the TestPyPI package are already copied
under `src/smog3/data/`, so normal installed SMOG3 runtime does not depend on
the legacy repository layout.
