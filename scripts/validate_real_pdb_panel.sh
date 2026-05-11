#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUT_DIR="real_pdb_validation"
DOCKER_IMAGE="smogserver/smog2:stable"
CASES="protein_ci2,protein_sh3,dna,rna,mixed_1a01_amp"
USE_INSTALLED_SMOG3=0
KEEP_EXISTING_DOWNLOADS=0
LOCAL_ONLY=0

usage() {
  cat <<'EOF'
Usage: bash scripts/validate_real_pdb_panel.sh [options]

Options:
  --use-installed-smog3       Run the installed smog3 console command instead of PYTHONPATH=src python3 -m smog3.smog2_native.
  --cases LIST                Comma-separated case names.
  --keep-existing-downloads   Reuse files already present in real_pdb_validation/inputs.
  --local-only                Skip RCSB downloads and use SMOG-CHECK fallback PDB files.
  --out-dir DIR               Output directory. Default: real_pdb_validation.
  -h, --help                  Show this help.

Default cases:
  protein_ci2,protein_sh3,dna,rna,mixed_1a01_amp

Additional cases:
  protein_ci2_ca,protein_sh3_ca
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --use-installed-smog3)
      USE_INSTALLED_SMOG3=1
      shift
      ;;
    --cases)
      CASES="${2:-}"
      shift 2
      ;;
    --keep-existing-downloads)
      KEEP_EXISTING_DOWNLOADS=1
      shift
      ;;
    --local-only)
      LOCAL_ONLY=1
      shift
      ;;
    --out-dir)
      OUT_DIR="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$CASES" ]]; then
  echo "--cases must not be empty" >&2
  exit 2
fi

case_spec() {
  case "$1" in
    protein_ci2) echo "protein_ci2|2CI2|2ci2_v2.pdb|-AA" ;;
    protein_ci2_ca) echo "protein_ci2_ca|2CI2|2ci2_v2.pdb|-CA" ;;
    protein_sh3) echo "protein_sh3|1SHG|1AKEapo_v2.pdb|-AA" ;;
    protein_sh3_ca) echo "protein_sh3_ca|1SHG|1AKEapo_v2.pdb|-CA" ;;
    dna) echo "dna|1BNA|DNA.terminal.pdb|-AA" ;;
    rna) echo "rna|1EHZ|tRNA.pdb|-AA" ;;
    mixed_1a01_amp) echo "mixed_1a01_amp|1A01|1A01-AMP.pdb|-AA" ;;
    *)
      echo "Unknown real-PDB panel case: $1" >&2
      return 1
      ;;
  esac
}

download_or_fallback() {
  local case_name="$1"
  local pdb_id="$2"
  local fallback_name="$3"
  local out_pdb="$4"
  local source_file="$5"
  local fallback="$ROOT/SMOG-CHECK/share/PDB.files/$fallback_name"

  if [[ "$KEEP_EXISTING_DOWNLOADS" -eq 1 && -s "$out_pdb" ]]; then
    echo "existing" > "$source_file"
    return 0
  fi

  if [[ "$LOCAL_ONLY" -eq 0 ]]; then
    local url="https://files.rcsb.org/download/${pdb_id}.pdb"
    local tmp="${out_pdb}.tmp"
    if command -v curl >/dev/null 2>&1 && curl -fsSL "$url" -o "$tmp"; then
      mv "$tmp" "$out_pdb"
      echo "download:${url}" > "$source_file"
      return 0
    fi
    rm -f "$tmp"
    echo "Download failed for $case_name ($url); falling back to $fallback_name" >&2
  fi

  if [[ ! -f "$fallback" ]]; then
    echo "Fallback PDB missing: $fallback" >&2
    return 1
  fi
  cp "$fallback" "$out_pdb"
  echo "fallback:$fallback_name" > "$source_file"
}

run_smog2_baseline() {
  local input_rel="$1"
  local out_rel="$2"
  local mode="$3"
  local log_path="$4"
  local -a mode_args
  read -r -a mode_args <<< "$mode"
  docker run --rm -v "$ROOT":/workdir "$DOCKER_IMAGE" bash -lc \
    'set -euo pipefail; cd /workdir; smog2 "$@"' _ \
    -i "/workdir/$input_rel" \
    "${mode_args[@]}" \
    -keep4SCM \
    -o "/workdir/$out_rel/model.top" \
    -g "/workdir/$out_rel/model.gro" \
    -n "/workdir/$out_rel/model.ndx" \
    -s "/workdir/$out_rel/model.contacts" \
    >"$log_path" 2>&1
}

run_smog3_candidate() {
  local input_rel="$1"
  local out_dir="$2"
  local mode="$3"
  local log_path="$4"
  local perl_log="$5"
  local no_perl_bin="$6"
  local -a mode_args
  read -r -a mode_args <<< "$mode"
  if [[ "$USE_INSTALLED_SMOG3" -eq 1 ]]; then
    env \
      PATH="$no_perl_bin:$PATH" \
      SMOG3_PERL_SENTINEL_LOG="$perl_log" \
      SMOG3_LEGACY_PERL_FALLBACK=0 \
      SMOG3_USE_SCM_DEFAULTS=1 \
      smog3 \
      -i "$input_rel" \
      "${mode_args[@]}" \
      -o "$out_dir/model.top" \
      -g "$out_dir/model.gro" \
      -n "$out_dir/model.ndx" \
      -s "$out_dir/model.contacts" \
      >"$log_path" 2>&1
  else
    env \
      PATH="$no_perl_bin:$PATH" \
      PYTHONPATH="$ROOT/src:${PYTHONPATH:-}" \
      SMOG3_PERL_SENTINEL_LOG="$perl_log" \
      SMOG3_LEGACY_PERL_FALLBACK=0 \
      SMOG3_USE_SCM_DEFAULTS=1 \
      python3 -m smog3.smog2_native \
      -i "$input_rel" \
      "${mode_args[@]}" \
      -o "$out_dir/model.top" \
      -g "$out_dir/model.gro" \
      -n "$out_dir/model.ndx" \
      -s "$out_dir/model.contacts" \
      >"$log_path" 2>&1
  fi
}

smog3_version() {
  if [[ "$USE_INSTALLED_SMOG3" -eq 1 ]]; then
    python3 -c 'import smog3; print(getattr(smog3, "__version__", "unknown"))' 2>/dev/null || echo "unknown"
  else
    PYTHONPATH="$ROOT/src" python3 -c 'import smog3; print(getattr(smog3, "__version__", "unknown"))' 2>/dev/null || echo "unknown"
  fi
}

smog2_docker_version() {
  local out
  out="$(docker run --rm "$DOCKER_IMAGE" bash -lc 'smog2 -v 2>&1 | head -20' 2>/dev/null || true)"
  if [[ "$out" =~ SMOG[[:space:]]v[0-9][^[:space:],]* ]]; then
    echo "${BASH_REMATCH[0]}"
  else
    echo "$out" | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g'
  fi
}

mkdir -p "$OUT_DIR"
if [[ "$KEEP_EXISTING_DOWNLOADS" -eq 0 ]]; then
  rm -rf "$OUT_DIR/inputs"
fi
rm -rf "$OUT_DIR/smog2_baseline" "$OUT_DIR/smog3_candidate" "$OUT_DIR/reports"
rm -f "$OUT_DIR/real_pdb_panel_summary.json" "$OUT_DIR/real_pdb_panel_summary.md"
mkdir -p "$OUT_DIR/inputs" "$OUT_DIR/smog2_baseline" "$OUT_DIR/smog3_candidate" "$OUT_DIR/reports"

NO_PERL_BIN="$OUT_DIR/no-perl-bin"
PERL_LOG="$OUT_DIR/smog3-perl-invocations.log"
mkdir -p "$NO_PERL_BIN"
cat > "$NO_PERL_BIN/perl" <<'EOF'
#!/usr/bin/env bash
echo "SMOG3 attempted to invoke perl: $*" >> "${SMOG3_PERL_SENTINEL_LOG:-/tmp/smog3-real-panel-perl.log}"
exit 127
EOF
chmod +x "$NO_PERL_BIN/perl"
: > "$PERL_LOG"

IFS=',' read -r -a requested_cases <<< "$CASES"
PANEL_RC=0

for raw_case in "${requested_cases[@]}"; do
  case_name="$(echo "$raw_case" | xargs)"
  [[ -n "$case_name" ]] || continue
  spec="$(case_spec "$case_name")" || exit 2
  IFS='|' read -r resolved_name pdb_id fallback_name mode <<< "$spec"
  input_pdb="$OUT_DIR/inputs/${resolved_name}.pdb"
  source_file="$OUT_DIR/inputs/${resolved_name}.source"
  baseline_dir="$OUT_DIR/smog2_baseline/$resolved_name"
  candidate_dir="$OUT_DIR/smog3_candidate/$resolved_name"
  report_json="$OUT_DIR/reports/${resolved_name}.json"
  baseline_log="$baseline_dir/smog2.log"
  candidate_log="$candidate_dir/smog3.log"
  mkdir -p "$baseline_dir" "$candidate_dir"

  echo "[real-panel] preparing $resolved_name ($pdb_id, $mode)"
  if ! download_or_fallback "$resolved_name" "$pdb_id" "$fallback_name" "$input_pdb" "$source_file"; then
    PYTHONPATH="$ROOT/src" python3 scripts/compare_real_pdb_panel.py compare \
      --case "$resolved_name" \
      --input-pdb "$input_pdb" \
      --input-source "download-error" \
      --baseline-dir "$baseline_dir" \
      --candidate-dir "$candidate_dir" \
      --baseline-rc 99 \
      --candidate-rc 99 \
      --baseline-log "$baseline_log" \
      --candidate-log "$candidate_log" \
      --report-json "$report_json" || true
    PANEL_RC=1
    continue
  fi

  input_source="$(cat "$source_file")"
  input_rel="${input_pdb#$ROOT/}"
  if [[ "$input_rel" == "$input_pdb" ]]; then
    input_rel="$input_pdb"
  fi

  echo "[real-panel] running official SMOG2 Docker baseline for $resolved_name"
  BASELINE_RC=0
  run_smog2_baseline "$input_rel" "${baseline_dir#$ROOT/}" "$mode" "$baseline_log" || BASELINE_RC=$?

  if [[ "$BASELINE_RC" -ne 0 && "$input_source" == download:* ]]; then
    echo "[real-panel] SMOG2 rejected downloaded $pdb_id for $resolved_name; retrying local fallback $fallback_name" >&2
    cp "$ROOT/SMOG-CHECK/share/PDB.files/$fallback_name" "$input_pdb"
    echo "fallback-after-smog2-error:$fallback_name" > "$source_file"
    input_source="$(cat "$source_file")"
    rm -rf "$baseline_dir"
    mkdir -p "$baseline_dir"
    baseline_log="$baseline_dir/smog2.log"
    BASELINE_RC=0
    run_smog2_baseline "$input_rel" "${baseline_dir#$ROOT/}" "$mode" "$baseline_log" || BASELINE_RC=$?
  fi

  CANDIDATE_RC=99
  if [[ "$BASELINE_RC" -eq 0 ]]; then
    echo "[real-panel] running SMOG3 candidate for $resolved_name"
    CANDIDATE_RC=0
    run_smog3_candidate "$input_rel" "$candidate_dir" "$mode" "$candidate_log" "$PERL_LOG" "$NO_PERL_BIN" || CANDIDATE_RC=$?
  else
    echo "[real-panel] skipping SMOG3 candidate for $resolved_name because SMOG2 baseline failed" >&2
    : > "$candidate_log"
  fi

  COMPARE_RC=0
  PYTHONPATH="$ROOT/src" python3 scripts/compare_real_pdb_panel.py compare \
    --case "$resolved_name" \
    --input-pdb "$input_pdb" \
    --input-source "$input_source" \
    --baseline-dir "$baseline_dir" \
    --candidate-dir "$candidate_dir" \
    --baseline-rc "$BASELINE_RC" \
    --candidate-rc "$CANDIDATE_RC" \
    --baseline-log "$baseline_log" \
    --candidate-log "$candidate_log" \
    --report-json "$report_json" || COMPARE_RC=$?
  if [[ "$COMPARE_RC" -ne 0 ]]; then
    PANEL_RC=1
  fi
done

PERL_COUNT=0
if [[ -s "$PERL_LOG" ]]; then
  PERL_COUNT="$(wc -l < "$PERL_LOG" | tr -d ' ')"
fi

SUMMARY_RC=0
PYTHONPATH="$ROOT/src" python3 scripts/compare_real_pdb_panel.py summarize \
  --reports-dir "$OUT_DIR/reports" \
  --summary-json "$OUT_DIR/real_pdb_panel_summary.json" \
  --summary-md "$OUT_DIR/real_pdb_panel_summary.md" \
  --smog3-version "$(smog3_version)" \
  --smog2-version "$(smog2_docker_version)" \
  --docker-image "$DOCKER_IMAGE" \
  --smog3-perl-invocations "$PERL_COUNT" \
  --smog3-perl-log "$PERL_LOG" || SUMMARY_RC=$?

echo "[real-panel] summary: $OUT_DIR/real_pdb_panel_summary.json"
cat "$OUT_DIR/real_pdb_panel_summary.md"

if [[ "$PERL_COUNT" -ne 0 ]]; then
  echo "[real-panel] ERROR: SMOG3 attempted to invoke Perl; see $PERL_LOG" >&2
  exit 3
fi

if [[ "$PANEL_RC" -ne 0 || "$SUMMARY_RC" -ne 0 ]]; then
  exit 1
fi
