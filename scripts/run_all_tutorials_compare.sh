#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Public tutorial comparison driver.  SMOG2 baselines are generated only inside
# smogserver/smog2:stable; SMOG3 candidates run from the source tree by default
# or from an installed `smog3` command with --use-installed-smog3.
PYTHONPATH="$ROOT/src:${PYTHONPATH:-}" python3 -m smog3.tutorial_validation "$@"
