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

Current official-image caveat: `smogserver/smog2:stable` reports SMOG v2.4.5
and cannot generate baselines for SMOG-CHECK cases 104 and 110-113 from this
repo's current testlist/templates. The image lacks the newer custom/pcos/ncos
dihedral function support needed by those templates, and `-warn -1` still exits
before complete baseline outputs are written. The parity runner leaves those
as `BASELINE_ERROR`; it does not mark them as passing or compare them against
candidate files.

You can also run comparator directly:

```bash
PYTHONPATH=src python3 -m smog3.parity_direct   --compare-existing   --baseline-dir parity_runs/case1/baseline   --candidate-dir parity_runs/case1/candidate   --report-json parity_case1.json
```
