#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

bash scripts/validate_real_pdb_panel.sh --local-only --cases protein_ci2
