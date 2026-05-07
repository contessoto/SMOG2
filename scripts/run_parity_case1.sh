#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

PYTHONPATH=src python -m smog3.parity_direct --cases 1 --report-json parity_case1.json
