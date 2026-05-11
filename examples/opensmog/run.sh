#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f 2CI2.pdb ]]; then
  if command -v curl >/dev/null 2>&1 && curl -fsSLO https://files.rcsb.org/download/2CI2.pdb; then
    :
  elif [[ -f ../../SMOG-CHECK/share/PDB.files/2ci2_v2.pdb ]]; then
    cp ../../SMOG-CHECK/share/PDB.files/2ci2_v2.pdb 2CI2.pdb
  else
    echo "Could not download 2CI2.pdb and no source-checkout fallback was found." >&2
    exit 1
  fi
fi

smog3 -i 2CI2.pdb -AA -OpenSMOG -OpenSMOGxml model.xml \
  -o model_os.top \
  -g model_os.gro \
  -n model_os.ndx \
  -s model_os.contacts

ls -lh model_os.top model_os.gro model_os.ndx model_os.contacts model.xml
