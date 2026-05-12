# Public Tutorial Validation

SMOG3 has a dedicated validation suite for the public SMOG tutorial/model categories listed at <https://smog-server.org/tutorials/>. This suite is separate from SMOG-CHECK: it is intended to make public-facing model examples reproducible against the official SMOG2 Docker image.

## Download Public Tutorial Assets

The repository does not commit downloaded tutorial inputs.  Fetch the public
tutorial pages and small linked assets with:

```bash
python3 scripts/fetch_smog_tutorial_assets.py
```

The fetcher starts at <https://smog-server.org/tutorials/>, saves HTML pages,
downloads small linked model-generation files, and writes an index:

```text
validation/tutorials/assets/
  pages/
  downloads/
  manifest_raw.json
```

Files larger than 100 MB are skipped as `LARGE_FILE_SKIPPED`.  Downloaded pages
and assets are ignored by git; the raw JSON index is kept so the crawl result can
be inspected.

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

Tutorial categories that require additional public tutorial input bundles are present in the manifest and classified as `MISSING_DOWNLOAD` until those files are fetched, or `MANUAL_INPUT_REQUIRED` when the public tutorial uses custom template/script steps that have not yet been mapped to an automated SMOG3 command.

## Commands

List manifest entries:

```bash
bash scripts/run_all_tutorials_compare.sh --list
```

Run implemented tutorial model-generation cases:

```bash
bash scripts/run_all_tutorials_compare.sh
```

Run every manifest entry, including `MISSING_INPUT` and `NOT_GENERATION_TEST` classifications:

```bash
bash scripts/run_all_tutorials_compare.sh --all
```

Run one case:

```bash
bash scripts/run_all_tutorials_compare.sh --case standard_aa_ci2
```

Fetch assets first, then run all manifest entries:

```bash
bash scripts/run_all_tutorials_compare.sh --all --download-first
```

Run only the official SMOG2 baseline side or only the SMOG3 side:

```bash
bash scripts/run_all_tutorials_compare.sh --case standard_aa_ci2 --smog2-only
bash scripts/run_all_tutorials_compare.sh --case standard_aa_ci2 --smog3-only
```

Run against an installed `smog3` command, such as the TestPyPI package:

```bash
PATH="/tmp/smog3-tutorial-test/bin:$PATH" \
bash scripts/run_all_tutorials_compare.sh --all --use-installed-smog3
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
bash scripts/run_all_tutorials_compare.sh --all --use-installed-smog3
```

The older `scripts/run_tutorial_validation.sh` wrapper is retained as a
compatibility alias for the same Python runner.

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
  tutorial_compare_summary.json
  tutorial_compare_summary.md
```

Per-case reports include command lines, return codes, file comparison results, topology section counts, contact counts, coordinate atom counts, XML parse status, and a first diff snippet when files differ.

## Status Meanings

- `PASS`: SMOG2 Docker and SMOG3 generated matching model-generation outputs
  under the documented comparator policy, or the requested one-sided generation
  check completed.
- `DIFF`: both sides ran, but at least one compared output differed.
- `SMOG2_ERROR`: the official Docker baseline command failed.
- `SMOG3_ERROR`: the SMOG3 candidate command failed.
- `MISSING_INPUT`: a required local checkout file is absent.
- `MISSING_DOWNLOAD`: a required public tutorial asset has not been fetched.
- `MANUAL_INPUT_REQUIRED`: public assets exist or are discoverable, but the
  tutorial includes custom preprocessing/template/script choices that still need
  explicit automated command mapping.
- `NOT_GENERATION_TEST`: the tutorial is a simulation/checkpoint/minimization
  workflow rather than a SMOG2 model-generation comparison.
- `UNSUPPORTED_BY_SMOG3`: the tutorial maps to a known unsupported SMOG3 feature.

## Comparator Policy

The tutorial suite reuses the same documented parity comparator used for SMOG-CHECK validation. It performs exact comparison except for harmless generated topology/XML metadata and tiny topology floating-point print noise already accepted by the SMOG-CHECK parity policy. It does not hide topology section value differences, contact differences, or coordinate differences.

## Perl Boundary

The SMOG2 baseline runs inside `smogserver/smog2:stable`, where legacy Perl is expected. SMOG3 candidate runs set `SMOG3_LEGACY_PERL_FALLBACK=0` and install a temporary `perl` sentinel first in `PATH`. Any SMOG3 runtime Perl attempt is logged and fails the tutorial validation run.
