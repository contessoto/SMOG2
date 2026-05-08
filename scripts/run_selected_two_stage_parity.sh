#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "pyproject.toml" ]]; then
  echo "Run from repository root" >&2
  exit 1
fi

PYTHONPATH=src python3 -m smog3.selected_parity \
  --cases "1,21,41,50,56,94" \
  --out-root "parity_runs/selected" \
  --report-json "parity_selected.json"
