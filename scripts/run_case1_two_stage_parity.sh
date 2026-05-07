#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "pyproject.toml" ]]; then
  echo "Run from repository root" >&2
  exit 1
fi

BASE_DIR="parity_runs/case1/baseline"
CAND_DIR="parity_runs/case1/candidate"

rm -rf "$BASE_DIR" "$CAND_DIR"
mkdir -p "$BASE_DIR" "$CAND_DIR"

echo "[1/3] Running baseline SMOG2 in official container"
docker run --rm -v "$PWD":/workdir smogserver/smog2:stable bash -lc '
set -euo pipefail
cd /workdir
mkdir -p parity_runs/case1/baseline
smog2 -i /workdir/SMOG-CHECK/share/PDB.files/1A01-AMP.pdb -AA -o /workdir/parity_runs/case1/baseline/model.top -g /workdir/parity_runs/case1/baseline/model.gro -n /workdir/parity_runs/case1/baseline/model.ndx -s /workdir/parity_runs/case1/baseline/model.contacts
'

echo "[2/3] Running candidate SMOG3 with local Python"
if ! command -v python3 >/dev/null 2>&1; then
  echo "Missing local python3 for candidate stage. Install python3 and retry." >&2
  exit 2
fi

CAND_CMD=(python3 -c "from smog3.smog2_native import main; import sys; raise SystemExit(main(sys.argv[1:]))" \
  -i SMOG-CHECK/share/PDB.files/1A01-AMP.pdb \
  -AA \
  -o "$CAND_DIR/model.top" \
  -g "$CAND_DIR/model.gro" \
  -n "$CAND_DIR/model.ndx" \
  -s "$CAND_DIR/model.contacts")

printf 'Candidate command: %q ' "${CAND_CMD[@]}"; echo
PYTHONPATH=src "${CAND_CMD[@]}"

MISSING=0
for f in model.top model.gro model.ndx model.contacts; do
  if [[ ! -f "$CAND_DIR/$f" ]]; then
    echo "ERROR: candidate output missing: $CAND_DIR/$f" >&2
    MISSING=1
  fi
done
if [[ "$MISSING" -ne 0 ]]; then
  echo "Candidate directory contents:" >&2
  ls -la "$CAND_DIR" >&2 || true
  find "$CAND_DIR" -maxdepth 2 -type f -print >&2 || true
  exit 3
fi

echo "[3/3] Comparing baseline vs candidate"
PYTHONPATH=src python3 -m smog3.parity_direct \
  --compare-existing \
  --baseline-dir "$BASE_DIR" \
  --candidate-dir "$CAND_DIR" \
  --report-json parity_case1.json

echo "Parity report written to parity_case1.json"
