# SMOG3

SMOG3 is a Python-native successor runtime for SMOG2 structure-based model
generation.  It keeps compatibility with the SMOG2 command behavior that is
covered by the repository validation suite while replacing normal SMOG3 runtime
execution with Python code.

SMOG3 is alpha software.  The current implementation has passed the original
drop-in SMOG-CHECK harness for tests 1 through 115 using a temporary `smog2`
wrapper that routes every SMOG-CHECK `smog2` call to Python SMOG3 code.

## Runtime Model

- SMOG3 normal runtime uses Python plus Java SCM where contact generation needs
  the existing `SCM.jar` tool.
- SMOG3 normal runtime does not invoke Perl.
- The legacy SMOG-CHECK harness itself is Perl.  That is acceptable for
  validation; the harness calls a temporary `smog2` wrapper that executes
  `python3 -m smog3.smogcheck_dropin_smog2`.
- Original SMOG2 baselines are generated only inside the official
  `smogserver/smog2:stable` Docker image.

## Install From TestPyPI

The next TestPyPI build prepared from this checkout is `0.1.0a2`. It includes
the tutorial validation suite and the expanded `smog3 --help` text that shows
OpenSMOG flags.

```bash
python3 -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  smog3
```

Verify the installed command:

```bash
smog3 --help
```

The help output lists the supported model flags and common generation options,
including `-OpenSMOG -OpenSMOGxml model.xml`.

## Basic Usage

Download a small PDB:

```bash
curl -fsSLO https://files.rcsb.org/download/2CI2.pdb
```

Generate an all-atom model:

```bash
smog3 -i 2CI2.pdb -AA \
  -o model.top \
  -g model.gro \
  -n model.ndx \
  -s model.contacts
```

Generate a C-alpha model:

```bash
smog3 -i 2CI2.pdb -CA \
  -o model_ca.top \
  -g model_ca.gro \
  -n model_ca.ndx \
  -s model_ca.contacts
```

Generate OpenSMOG XML:

```bash
smog3 -i 2CI2.pdb -AA -OpenSMOG -OpenSMOGxml model.xml \
  -o model_os.top \
  -g model_os.gro \
  -n model_os.ndx \
  -s model_os.contacts
```

## Validation

From a source checkout:

```bash
PYTHONPATH=src pytest -q tests
bash scripts/run_selected_two_stage_parity.sh
bash scripts/run_smogcheck_dropin_smog3.sh 1 115
```

Or run the combined source-checkout validation script:

```bash
bash scripts/run_all_tests.sh
```

Real-world panel comparison against official SMOG2 Docker:

```bash
bash scripts/validate_real_pdb_panel.sh
bash scripts/validate_real_pdb_panel.sh --use-installed-smog3
```

Public SMOG tutorial/model-generation validation:

```bash
bash scripts/run_tutorial_validation.sh --list
bash scripts/run_tutorial_validation.sh --all
bash scripts/run_tutorial_validation.sh --case standard_aa_ci2
bash scripts/run_tutorial_validation.sh --all --use-installed-smog3
```

The tutorial suite compares SMOG3-generated output files against official SMOG2
Docker baselines for representative public tutorial categories. It reports
`PASS` only when `.top`, `.gro`/`.g96`, `.ndx`, `.contacts`, and `.xml` outputs
match under the documented parity comparator policy.

The validation outputs are written to ignored local directories such as
`parity_runs/`, `smogcheck_dropin_runs/`, `real_pdb_validation/`, and
`validation/tutorials/runs/`.

## Repository Layout

- `src/smog3/`: Python-native SMOG3 package.
- `src/smog3/data/`: packaged runtime templates and `SCM.jar`.
- `tests/`: Python tests for native tools, parity behavior, and validation
  helpers.
- `scripts/`: validation and release helper scripts.
- `docs/`: installation, quickstart, validation, TestPyPI, and developer notes.
- `examples/`: source-friendly examples that do not commit generated outputs.
- `validation/`: notes for SMOG-CHECK and real PDB validation campaigns.
- `legacy/smog2/`: notes about the preserved legacy SMOG2 layout.
- `SMOG-CHECK/`, `share/`, and legacy files at the repository root are retained
  because the proven validation harnesses still depend on those historical
  paths.

## License And Citation

SMOG3 is distributed under the same GPL-compatible repository license as SMOG2.
See `COPYING` and `CITATION.cff`.
