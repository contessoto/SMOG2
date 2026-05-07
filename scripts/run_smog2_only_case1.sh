#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "pyproject.toml" ]]; then
  echo "Run this script from repository root (pyproject.toml not found)." >&2
  exit 1
fi

mkdir -p baseline_case1

TTY_FLAGS=""
if [[ -t 0 && -t 1 ]]; then
  TTY_FLAGS="-it"
fi

docker run --rm ${TTY_FLAGS} -v "$PWD":/workdir smogserver/smog2:stable bash -lc '
set -euo pipefail
cd /workdir
which smog2
smog2 -v

mkdir -p baseline_case1
cd baseline_case1
smog2 -i /workdir/SMOG-CHECK/share/PDB.files/1A01-AMP.pdb -AA -o model.top -g model.gro -n model.ndx -s model.contacts
ls -1 model.top model.gro model.ndx model.contacts
'

echo "Baseline outputs written under ./baseline_case1"
