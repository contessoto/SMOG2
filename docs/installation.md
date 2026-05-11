# Installation

## Install The Alpha From TestPyPI

```bash
python3 -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  smog3
```

Verify the command:

```bash
smog3 --help
python3 -c "import smog3; print(smog3.__version__)"
```

## Editable Source Install

```bash
git clone https://github.com/smog-server/SMOG2.git
cd SMOG2
python3 -m pip install -e .
PYTHONPATH=src pytest -q tests
```

## Runtime Requirements

- Python 3.9 or newer.
- Java on `PATH` for SCM contact generation.  SMOG3 bundles the `SCM.jar`
  runtime file in the Python package data.
- Perl is not required for normal SMOG3 runtime.

## Optional Comparison Requirements

Docker is required only when comparing with original SMOG2 through the official
baseline image:

```bash
docker pull smogserver/smog2:stable
```

The validation scripts use Docker to run original SMOG2 and then run SMOG3
locally.

## Troubleshooting

- If contact files are missing, check that `java` is available:
  `java -version`.
- If Docker baseline runs fail, check `docker info` and make sure
  `smogserver/smog2:stable` can be pulled.
- If an installed `smog3` command is not found, check the active virtual
  environment and `python3 -m pip show smog3`.
- If the original SMOG-CHECK harness prints `SMOG2`, remember that the harness
  text is legacy wording; the drop-in wrapper logs each call routed to SMOG3.
