# SMOG2 vs SMOG3 Direct Parity Environment

This project includes a reproducible local container/devcontainer setup for running baseline SMOG2 (Perl) and SMOG3 (Python-native) parity checks.

## Included baseline dependencies

The devcontainer installs:

- `perl`
- `libxml-simple-perl`
- `libxml-validator-schema-perl`
- `pdl`
- Java runtime (`default-jre-headless`) for SCM/WHAM tooling needs
- Python 3 + pip/build tools

It also installs SMOG3 in editable mode.

## Build and run with Docker

From repository root:

```bash
docker build -f .devcontainer/Dockerfile -t smog3-parity .
docker run --rm -it -v "$PWD":/workspace/SMOG2 smog3-parity bash
```

Inside container, verify baseline Perl deps:

```bash
perl -MXML::Simple -MXML::Validator::Schema -MPDL -e 'print "SMOG2 Perl deps ok\n"'
```

## Run case-1 direct parity harness

```bash
./scripts/run_parity_case1.sh
```

This runs:

```bash
PYTHONPATH=src python -m smog3.parity_direct --cases 1 --report-json parity_case1.json
```

## Optional: run test suite

```bash
PYTHONPATH=src pytest -q tests
```

## VS Code devcontainer

Open this repo in VS Code and choose **Reopen in Container**. The configuration is in:

- `.devcontainer/devcontainer.json`
- `.devcontainer/Dockerfile`


## GitHub Actions parity

When local Docker is unavailable, use the CI workflow:

- `.github/workflows/smog3-parity.yml`

It installs Perl deps (`XML::Simple`, `XML::Validator::Schema`, `PDL`), runs `pytest`, runs direct parity case 1, and uploads `parity_case1.json` as an artifact even on failure.


## Official SMOG2 Docker baseline (recommended on macOS)

The official SMOG2 docs and image:

- Docs: https://smog-server.org/smog2/docker/
- Image: `smogserver/smog2:stable`

### 1) Install Docker Desktop

Install Docker Desktop for macOS and ensure `docker` works in Terminal.

### 2) Pull official image

```bash
docker pull smogserver/smog2:stable
```

### 3) Optional sanity-check shell

```bash
docker run -it --rm -v "$PWD":/workdir smogserver/smog2:stable
```

Inside container, `/workdir` is your mounted repo.

### 4) Run baseline-only case 1 (SMOG2 only)

```bash
./scripts/run_smog2_only_case1.sh
```

This writes baseline outputs to `./baseline_case1`.

### 5) Run full official-container parity (SMOG2 baseline + SMOG3 candidate)

```bash
./scripts/run_official_smog2_docker_parity.sh
```

This performs:
- `which smog2`
- `smog2 -v`
- `test -d /opt/smog2/SMOG-CHECK`
- Python/pip checks
- editable install of SMOG3
- `PYTHONPATH=src python3 -m pytest -q tests`
- `PYTHONPATH=src python3 -m smog3.parity_direct --cases 1 --report-json parity_case1.json`

If parity fails, send `parity_case1.json` back to Codex for targeted diff fixes.


## Two-stage parity workflow (official SMOG2 Docker + local Python SMOG3)

Because `smogserver/smog2:stable` may not include Python, use a two-stage flow:

1. Baseline stage inside official SMOG2 Docker image.
2. Candidate stage with local `python3` running SMOG3.
3. Compare existing outputs with `smog3.parity_direct --compare-existing`.

Run from repo root:

```bash
bash scripts/run_case1_two_stage_parity.sh
```

This script will:
- recreate `parity_runs/case1/baseline`
- recreate `parity_runs/case1/candidate`
- run baseline SMOG2 in `smogserver/smog2:stable`
- run candidate SMOG3 locally
- compare outputs and write `parity_case1.json`
- exit nonzero if parity fails

For `model.top`, the comparator ignores only the free-form header metadata
before the first topology section. Section contents and ordering remain
strict. Inside numeric topology parameter sections (`[ atomtypes ]`,
`[ bonds ]`, `[ angles ]`, `[ dihedrals ]`, and `[ pairs ]`), the comparator
also accepts tiny floating-point print ULP differences when row identity,
token layout, comments, and all non-floating tokens match exactly. The same
narrow rule treats `+180` and `-180` as equivalent only for the `[ dihedrals ]`
`phi0` column because they are the same periodic endpoint printed with
different signs. Larger numeric differences, reordered rows, changed atom
indices, changed function types, and changed comments remain failures.

For `model.xml`, the comparator similarly ignores only the generated XML
comment block before the `<OpenSMOGforces>` root element. All OpenSMOG force
elements, interaction ordering, and parameter values remain strict.

## Current SMOG-CHECK parity status

Latest verified campaign:

| Result | Count |
| --- | ---: |
| Total direct SMOG-CHECK cases | 115 |
| PASS | 110 |
| FAIL | 0 |
| UNSUPPORTED | 0 |
| BASELINE_ERROR | 5 |
| CANDIDATE_ERROR | 0 |

SMOG3 normal runtime uses zero Perl. Original SMOG2 is run only inside the
official `smogserver/smog2:stable` Docker image to generate baseline outputs.
SMOG3 uses Python plus Java where needed for SCM contact generation.

## Drop-in SMOG-CHECK harness compatibility

The two-stage parity campaign above is not the same as running the original
`SMOG-CHECK/smog-check` script with `smog2` replaced on `PATH`.  The original
harness has SMOG2 installation-layout assumptions: it reads the last line of
the `smog2` executable to infer `SMOG_PATH`, Perl module paths, template
directories, and `src/tools/SCM.jar`.

`scripts/run_smogcheck_dropin_smog3.sh` creates a temporary compatibility
layout under `smogcheck_dropin_runs/`, creates a temporary `bin/smog2` wrapper,
puts that directory first in `PATH`, and runs the original `SMOG-CHECK`
harness from a copied work directory.  The wrapper body invokes:

```bash
python3 -m smog3.smogcheck_dropin_smog2
```

It sets `SMOG3_LEGACY_PERL_FALLBACK=0` and logs every wrapper invocation to
`smog3-wrapper-invocations.jsonl`.  The final trailer line in the wrapper is
present only so the original harness can discover its expected SMOG2-style
layout; execution exits before that trailer and does not call Perl from SMOG3.
The harness itself is still Perl because `SMOG-CHECK/smog-check` is Perl-based.

Latest drop-in harness checks:

| Harness command | Result |
| --- | --- |
| `bash scripts/run_smogcheck_dropin_smog3.sh 1 50` | 50/50 PASS |
| `bash scripts/run_smogcheck_dropin_smog3.sh 51 60` | 10/10 PASS |
| `bash scripts/run_smogcheck_dropin_smog3.sh 61 80` | 20/20 PASS |
| `bash scripts/run_smogcheck_dropin_smog3.sh 81 115` | PASS through 93, OpenSMOG/DIHE failures afterward, harness crash at 112 |
| `bash scripts/run_smogcheck_dropin_smog3.sh 114 115` | 2/2 PASS |

Aggregating the direct harness runs above gives 99 PASS, 14 FAIL, and 2
HARNESS_CRASH cases.  This is intentionally reported separately from the
two-stage Docker-baseline parity result.

Observed direct-harness failures:

| Cases | Classification | Evidence |
| --- | --- | --- |
| 94-104, 106 | OpenSMOG output differences versus local SMOG-CHECK v2.7beta checks | The harness reports `TOP FIELDS FOUND` and `OPENSMOG: EXCLUSIONS`; for case 94 it expected XML exclusion generation `1` but found `0`. |
| 110-111 | DIHE template topology/defaults mismatch | The harness reports `DEFAULTS, num entries` because it expected 3 default fields and found 5. |
| 112-113 | SMOG-CHECK harness internal crash | The harness exits with `Internal error : DUPLICATE DIHEDRAL CHECKING` before completing the cases. |

These direct-harness failures are not hidden by the comparator and should not
be described as full SMOG-CHECK success.  They are the next compatibility
targets if the goal shifts from Docker-baseline output parity to satisfying the
local original SMOG-CHECK v2.7beta introspection checks directly.

The five `BASELINE_ERROR` cases are not counted as SMOG3 failures. In each
case, the official Docker image reports SMOG v2.4.5 and exits before complete
baseline outputs are written. `-warn -1` was also checked and still did not
produce complete baseline files.

| Case | SMOG-CHECK entry | Baseline command options | Official-image baseline failure | Later validation needs |
| ---: | --- | --- | --- | --- |
| 104 | `1F4N_v2 OpenSMOG AA-CCD default` | `-OpenSMOG -OpenSMOGxml model.xml -t SMOG-CHECK/share/templates/SBM_AA+customContacts+customDihedrals` | SIF schema validation rejects `OpenSMOGtype="dihedral"` before output files are complete. | A SMOG2 baseline image/release that supports the current custom dihedral OpenSMOG template schema. |
| 110 | `RNA+protein AA-DIHE default` | `-t SMOG-CHECK/share/templates/SBM_AA_DIHE` | SMOG v2.4.5 rejects `dihedral_ncos`/`dihedral_pcos` template functions as unsupported. | A SMOG2 baseline image/release with pcos/ncos dihedral support. |
| 111 | `RNA+protein OpenSMOG AA-DIHE default` | `-OpenSMOG -OpenSMOGxml model.xml -t SMOG-CHECK/share/templates/SBM_AA_DIHE` | SMOG v2.4.5 rejects `dihedral_ncos`/`dihedral_pcos` before complete topology/XML output. | A SMOG2 baseline image/release with pcos/ncos dihedral and OpenSMOG support for these templates. |
| 112 | `RNA+protein AA-DIHE4 default` | `-t SMOG-CHECK/share/templates/SBM_AA_DIHE4` | SMOG v2.4.5 rejects `dihedral_pcos4`/`dihedral_ncos4` template functions as unsupported. | A SMOG2 baseline image/release with pcos4/ncos4 dihedral support. |
| 113 | `RNA+protein OpenSMOG AA-DIHE4 default` | `-OpenSMOG -OpenSMOGxml model.xml -t SMOG-CHECK/share/templates/SBM_AA_DIHE4` | SMOG v2.4.5 rejects `dihedral_ncos4`/`dihedral_pcos4` before complete topology/XML output. | A SMOG2 baseline image/release with pcos4/ncos4 dihedral and OpenSMOG support for these templates. |

## User commands

Run the Python test suite:

```bash
PYTHONPATH=src pytest -q tests
```

Run the selected representative parity set:

```bash
bash scripts/run_selected_two_stage_parity.sh
```

Run the original SMOG-CHECK harness against a SMOG3-backed `smog2` wrapper:

```bash
bash scripts/run_smogcheck_dropin_smog3.sh 1 50
```

Run the full SMOG-CHECK-style parity campaign:

```bash
bash scripts/run_all_smogcheck_two_stage_parity.sh --all
```

Inspect generated reports:

```bash
cat parity_selected.json
cat parity_all_summary.json
ls parity_runs/all/reports
```

You can also run comparator directly:

```bash
PYTHONPATH=src python3 -m smog3.parity_direct   --compare-existing   --baseline-dir parity_runs/case1/baseline   --candidate-dir parity_runs/case1/candidate   --report-json parity_case1.json
```
