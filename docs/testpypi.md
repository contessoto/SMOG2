# TestPyPI Release

## Build

```bash
python3 -m pip install -U build twine
noglob rm -rf dist build *.egg-info src/*.egg-info
python3 -m build
```

## Check

```bash
python3 -m twine check dist/*
```

## Upload To TestPyPI

Manual upload:

```bash
python3 -m twine upload --repository testpypi dist/*
```

Preferred trusted-publishing path:

1. Configure TestPyPI trusted publishing for this repository.
2. Use the workflow `.github/workflows/publish-testpypi.yml`.
3. Run it manually with `workflow_dispatch`.

The workflow builds distributions, runs `twine check`, and publishes to:

```text
https://test.pypi.org/legacy/
```

## Install From TestPyPI

```bash
python3 -m venv /tmp/smog3-testpypi
/tmp/smog3-testpypi/bin/python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  smog3
```

## Verify

```bash
/tmp/smog3-testpypi/bin/python -c "import smog3; print(smog3.__version__)"
/tmp/smog3-testpypi/bin/smog3 --help
```

Run a real-PDB comparison with the installed command:

```bash
PATH="/tmp/smog3-testpypi/bin:$PATH" \
  bash scripts/validate_real_pdb_panel.sh --use-installed-smog3 --local-only --cases protein_ci2
```
