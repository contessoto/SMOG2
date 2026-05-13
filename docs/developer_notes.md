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

## Runtime vs Validation-Only Code

Runtime code is the path reachable from the public installed commands:

- `smog3`: native model generation through `smog3.smog2_native.main`.
- `smog3-adjustpdb`: native PDB preprocessing for tutorial workflows.
- `smog3-extract`: native fragment extraction and index remapping.
- `smog3-ions`: native explicit-ion topology/GRO/XML post-processing.
- `smog3-parity-direct`: a validation helper, not required for normal model
  generation.

The normal `smog3` runtime must not call Perl, original SMOG2, Docker, or the
legacy `src/smogv2` script.  Validation-only modules and scripts may mention
Docker or SMOG2 because they generate official baselines:

- `smogcheck_parity.py`
- `tutorial_validation.py`
- `parity_direct.py`
- `scripts/run_*parity*.sh`
- `scripts/run_smogcheck_dropin_smog3.sh`

The drop-in SMOG-CHECK wrapper is a compatibility path for the legacy harness.
It sets `SMOG3_LEGACY_PERL_FALLBACK=0` before dispatching to native Python.

## Java SCM

SMOG3 uses Java for the existing SCM contact-generation workflow.  The Python
package bundles `SCM.jar` in `smog3/data/tools/SCM.jar`.

Java is optional only for workflows that do not require SCM contacts.  For
normal parity-covered model generation, install Java and verify it with:

```bash
java -version
```

## Package Data

Runtime templates and tools are packaged under `src/smog3/data/` and included
by `pyproject.toml`:

- `data/share/templates/`: default AA, CA, and Gaussian template files.
- `data/SMOG-CHECK/share/templates/`: validation-covered custom templates used
  to reproduce SMOG-CHECK behavior from an installed wheel.
- `data/tools/SCM.jar`: Java SCM contact generation.

When adding data needed by installed-wheel runtime, put it under
`src/smog3/data/` and verify a wheel install outside the repository can still
run AA, CA, and OpenSMOG examples.

## Adding a Model Flag

1. Add the command-line flag to `smog2_native.main`.
2. Map the flag to the appropriate template/model branch.
3. Implement topology/contact behavior in Python, using SMOG2 source/templates
   only as behavioral reference.
4. Add unit coverage for the new branch.
5. Add or extend two-stage parity coverage against official SMOG2 Docker.
6. Confirm the drop-in SMOG-CHECK wrapper translation if the legacy harness uses
   the new flag.

Do not call Perl or original SMOG2 from the runtime implementation.

## Adding a Tutorial Validation Case

1. Fetch public tutorial assets with `python3 scripts/fetch_smog_tutorial_assets.py`.
2. Add a case to `validation/tutorials/tutorial_manifest.yml`.
3. Include the exact public preprocessing/model-generation commands where
   possible, including `smog_adjustPDB`, `smog_ions`, maps, templates, and
   tutorial scripts.
4. Run `bash scripts/run_all_tutorials_compare.sh --case <case_id>`.
5. Keep the case as `NOT_GENERATION_TEST` only when it is simulation-only.

The tutorial runner compares official SMOG2 Docker outputs with SMOG3 outputs;
it should never mark a model-generation workflow PASS without generated files
being compared.

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

Run it with:

```bash
bash scripts/run_smogcheck_dropin_smog3.sh 1 115
```

The wrapper invocation log records each routed SMOG3 call.  Perl use by the
SMOG-CHECK harness itself is expected; Perl use inside SMOG3 runtime is not.

## TestPyPI Release Checklist

1. Update `pyproject.toml`, `src/smog3/__init__.py`, README, and CHANGELOG.
2. Clean old build outputs: `rm -rf dist build *.egg-info src/*.egg-info`.
3. Build: `python3 -m build`.
4. Check: `python3 -m twine check dist/*`.
5. Install the wheel in a fresh virtual environment.
6. Run AA, CA, and OpenSMOG generation from an empty directory.
7. Run `PYTHONPATH=src pytest -q tests`.
8. Upload only to TestPyPI until a real PyPI release is approved:

```bash
python3 -m twine upload --repository testpypi dist/*
```

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
