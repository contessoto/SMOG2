#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "pyproject.toml" ]]; then
  echo "Run from repository root" >&2
  exit 1
fi

PYTHONPATH=src python3 -m smog3.smogcheck_parity "$@"
