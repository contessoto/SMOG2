# Public Tutorial Validation

SMOG3 has a dedicated validation suite for the public SMOG tutorial/model categories listed at <https://smog-server.org/tutorials/>. This suite is separate from SMOG-CHECK: it is intended to make public-facing model examples reproducible against the official SMOG2 Docker image.

## What Is Tested

The tutorial runner validates SMOG model-generation outputs:

- `model.top`
- `model.gro` or `model.g96`
- `model.ndx`
- `model.contacts`
- `model.xml` for OpenSMOG cases

Simulation execution is not part of this suite. Tutorials that are only about OpenMM, Gromacs, checkpoint, or minimization execution are classified as `NOT_GENERATION_TEST`.

## Current Implemented Panel

The initial implemented tutorial panel uses local SMOG-CHECK inputs that correspond to public tutorial categories:

- standard all-atom CI2
- standard C-alpha CI2
- OpenSMOG all-atom CI2
- OpenSMOG C-alpha CI2
- all-atom ligand model with `1A01-AMP`
- all-atom DNA terminal model
- all-atom tRNA model
- all-atom disulfide local model
- C-alpha disulfide local model
- all-atom user-contact/multiple-contact local model

Tutorial categories that require additional public tutorial input bundles are present in the manifest and classified as `MISSING_INPUT` until those files are added.

## Commands

List manifest entries:

```bash
bash scripts/run_tutorial_validation.sh --list
```

Run implemented tutorial model-generation cases:

```bash
bash scripts/run_tutorial_validation.sh
```

Run every manifest entry, including `MISSING_INPUT` and `NOT_GENERATION_TEST` classifications:

```bash
bash scripts/run_tutorial_validation.sh --all
```

Run one case:

```bash
bash scripts/run_tutorial_validation.sh --case standard_aa_ci2
```

Run against an installed `smog3` command, such as the TestPyPI package:

```bash
PATH="/tmp/smog3-tutorial-test/bin:$PATH" \
bash scripts/run_tutorial_validation.sh --all --use-installed-smog3
```

## Installed TestPyPI Workflow

```bash
python3 -m venv /tmp/smog3-tutorial-test
/tmp/smog3-tutorial-test/bin/python -m pip install -U pip
/tmp/smog3-tutorial-test/bin/python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  smog3

PATH="/tmp/smog3-tutorial-test/bin:$PATH" \
bash scripts/run_tutorial_validation.sh --all --use-installed-smog3
```

## Reports

Each run writes:

```text
validation/tutorials/runs/<timestamp>/
  inputs/
  smog2_baseline/<case_id>/
  smog3_candidate/<case_id>/
  reports/<case_id>.json
  tutorial_validation_summary.json
  tutorial_validation_summary.md
```

Per-case reports include command lines, return codes, file comparison results, topology section counts, contact counts, coordinate atom counts, XML parse status, and a first diff snippet when files differ.

## Comparator Policy

The tutorial suite reuses the same documented parity comparator used for SMOG-CHECK validation. It performs exact comparison except for harmless generated topology/XML metadata and tiny topology floating-point print noise already accepted by the SMOG-CHECK parity policy. It does not hide topology section value differences, contact differences, or coordinate differences.

## Perl Boundary

The SMOG2 baseline runs inside `smogserver/smog2:stable`, where legacy Perl is expected. SMOG3 candidate runs set `SMOG3_LEGACY_PERL_FALLBACK=0` and install a temporary `perl` sentinel first in `PATH`. Any SMOG3 runtime Perl attempt is logged and fails the tutorial validation run.
