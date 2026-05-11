#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Thin wrapper around the Python tutorial validation runner. The runner keeps
# SMOG2 inside the official Docker image and routes SMOG3 through either the
# source tree or an installed smog3 console command.
PYTHONPATH="$ROOT/src:${PYTHONPATH:-}" python3 -m smog3.tutorial_validation "$@"
