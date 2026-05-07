#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "pyproject.toml" ]]; then
  echo "Run this script from repository root (pyproject.toml not found)." >&2
  exit 1
fi

TTY_FLAGS=""
if [[ -t 0 && -t 1 ]]; then
  TTY_FLAGS="-it"
fi

docker run --rm ${TTY_FLAGS} -v "$PWD":/workdir smogserver/smog2:stable bash -lc '
set -euo pipefail
cd /workdir

echo "[check] smog2 availability"
which smog2
smog2 -v

echo "[check] SMOG-CHECK path"
test -d /opt/smog2/SMOG-CHECK

echo "[check] Python availability"
if ! command -v python3 >/dev/null 2>&1; then
  echo "Missing python3 in smogserver/smog2:stable container." >&2
  exit 2
fi

if ! python3 -m pip --version >/dev/null 2>&1; then
  echo "Missing pip for python3 in smogserver/smog2:stable container." >&2
  echo "Install pip in the container (or use a derived image) before running parity." >&2
  exit 2
fi

python3 -m pip install -U pip
python3 -m pip install -e .
python3 -m pip install pytest

PYTHONPATH=src python3 -m pytest -q tests
PYTHONPATH=src python3 -m smog3.parity_direct --cases 1 --report-json parity_case1.json

echo "Parity report written to /workdir/parity_case1.json"
'
