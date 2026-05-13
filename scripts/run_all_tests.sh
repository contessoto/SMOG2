#!/usr/bin/env bash
set -euo pipefail

# Run the release-level SMOG3 validation checks from a source checkout.
# This intentionally exercises both the Python test suite and the legacy
# SMOG-CHECK harness routed through the Python-native SMOG3 drop-in wrapper.
if [[ ! -f "pyproject.toml" ]]; then
  echo "Run from repository root" >&2
  exit 1
fi

# Print a readable section header, run the command, and stop on the first
# failure so the command can be used directly in CI or release checks.
run_step() {
  local name="$1"
  shift
  echo
  echo "==> $name"
  "$@"
  echo "PASS: $name"
}

run_step "Python tests" env PYTHONPATH=src pytest -q tests
run_step "Selected two-stage parity" bash scripts/run_selected_two_stage_parity.sh
run_step "Drop-in SMOG-CHECK 1-115" bash scripts/run_smogcheck_dropin_smog3.sh 1 115
run_step "Public tutorial validation" bash scripts/run_all_tutorials_compare.sh --all --download-first

echo
echo "All SMOG3 validation steps passed."
