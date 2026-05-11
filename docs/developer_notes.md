# Developer Notes

## Architecture

SMOG3 keeps the Python package in `src/smog3`.  The main runtime entry point is
`smog3.smog2_native.main`, which implements the SMOG2-compatible generation
path used by the `smog3` console command and the drop-in SMOG-CHECK wrapper.

The implementation is intentionally parity-driven:

- original SMOG2 baselines are generated only inside
  `smogserver/smog2:stable`;
- SMOG3 candidates are generated locally with Python plus Java SCM;
- reports compare generated outputs under a documented metadata/tiny-print
  tolerance policy.

## `smog2_native.py`

`smog2_native.py` contains the current SMOG2-compatible native generator.  It
handles PDB parsing, all-atom and C-alpha model generation, topology writing,
GRO/NDX writing, OpenSMOG XML output, contact generation setup, and the case
families covered by SMOG-CHECK.

Runtime data is resolved first from packaged files under `src/smog3/data`, then
from the source checkout fallback layout.  That lets TestPyPI-installed SMOG3
run without depending on the repository root.

## Java SCM

SMOG3 uses Java for the existing SCM contact-generation workflow.  The Python
package bundles `SCM.jar` in `smog3/data/tools/SCM.jar`.

Java is optional only for workflows that do not require SCM contacts.  For
normal parity-covered model generation, install Java and verify it with:

```bash
java -version
```

## Parity Harnesses

- `scripts/run_selected_two_stage_parity.sh` runs representative two-stage
  parity cases.
- `scripts/run_all_smogcheck_two_stage_parity.sh` runs broader two-stage
  campaigns against the official Docker baseline.
- `scripts/run_smogcheck_dropin_smog3.sh` runs the original SMOG-CHECK harness
  with a temporary SMOG3-backed `smog2` wrapper.
- `scripts/validate_real_pdb_panel.sh` compares real or fallback PDBs against
  official SMOG2 Docker.

## Drop-In Wrapper

The drop-in wrapper exists because SMOG-CHECK expects an executable named
`smog2` and inspects wrapper text.  The wrapper has legacy-looking metadata
for the harness, but the executable path calls Python SMOG3 and sets
`SMOG3_LEGACY_PERL_FALLBACK=0`.

## Legacy Layout

The original SMOG2 files are preserved in their historical locations for now
because the proven validation harnesses depend on those paths.  The
`legacy/smog2` directory documents that decision.  A future migration can move
legacy assets after scripts and tests are adjusted together.

## Known Limitations

- SMOG3 is still alpha software.
- TestPyPI packages do not install a public `smog2` command by default.
- Full comparison against original SMOG2 requires Docker.
- The original SMOG-CHECK harness is Perl, even though SMOG3 runtime is not.
