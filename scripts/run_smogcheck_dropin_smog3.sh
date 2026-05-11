#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON:-python3}"
RUN_ROOT="${SMOG3_DROPIN_RUN_ROOT:-$ROOT/smogcheck_dropin_runs}"
STAMP="$(date +%Y%m%d-%H%M%S)"
mkdir -p "$RUN_ROOT"
RUN_DIR="$(mktemp -d "$RUN_ROOT/$STAMP.XXXXXX")"
WORK_DIR="$RUN_DIR/SMOG-CHECK"
BIN_DIR="$RUN_DIR/bin"
LAYOUT_DIR="$RUN_DIR/smog2-layout"
LOG_FILE="$RUN_DIR/smog3-wrapper-invocations.jsonl"
OUT_FILE="$RUN_DIR/smog-check.out"
PERL_BIN="$(command -v perl)"

mkdir -p "$BIN_DIR" "$LAYOUT_DIR"
cp -R "$ROOT/SMOG-CHECK" "$WORK_DIR"

ln -s "$ROOT/share" "$LAYOUT_DIR/share"
ln -s "$ROOT/src" "$LAYOUT_DIR/src"
for template in \
  "SBM_AA" \
  "SBM_AA+gaussian" \
  "SBM_calpha" \
  "SBM_calpha+gaussian" \
  "SBM_2cg" \
  "SBM_cr2" \
  "SBM_AA+customContacts" \
  "SBM_AA+customContacts+customDihedrals" \
  "SBM_AA_BOND" \
  "SBM_AA_DIHE" \
  "SBM_AA_DIHE4" \
  "SBM_match"
do
  ln -s "$ROOT/share/templates/$template" "$LAYOUT_DIR/$template"
done

cat > "$BIN_DIR/smog2" <<WRAPPER
#!/usr/bin/env bash
set -euo pipefail
export PERLLIB=$LAYOUT_DIR/src
export PERL5LIB=$LAYOUT_DIR/src
export PYTHONPATH="$ROOT/src:\${PYTHONPATH:-}"
export SMOG3_LEGACY_PERL_FALLBACK=0
export SMOG3_USE_SCM_DEFAULTS=1
export SMOG3_DROPIN_WRITE_USER_CONTACTS=1
export SMOG3_DROPIN_LOG="$LOG_FILE"
"$PYTHON_BIN" -m smog3.smogcheck_dropin_smog2 "\$@"
rc=\$?
exit "\$rc"
$PERL_BIN $LAYOUT_DIR/src/smogv2
WRAPPER
chmod +x "$BIN_DIR/smog2"

echo "Drop-in SMOG-CHECK run directory: $RUN_DIR"
echo "Drop-in smog2 wrapper: $BIN_DIR/smog2"
echo "SMOG2 layout advertised to SMOG-CHECK: $LAYOUT_DIR"
echo "Arguments passed to SMOG-CHECK: ${*:-<all>}"

set +e
(
  cd "$WORK_DIR"
  PATH="$BIN_DIR:$PATH" ./smog-check "$@"
) > "$OUT_FILE" 2>&1
rc=$?
set -e

cat "$OUT_FILE"

if [[ -s "$LOG_FILE" ]]; then
  invocations="$(wc -l < "$LOG_FILE" | tr -d ' ')"
else
  invocations="0"
fi
echo
echo "SMOG3 wrapper invocations: $invocations"
echo "Wrapper invocation log: $LOG_FILE"
echo "SMOG-CHECK transcript: $OUT_FILE"

exit "$rc"
