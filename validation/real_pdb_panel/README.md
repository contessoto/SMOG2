# Real PDB Panel Validation

The real-PDB panel downloads selected structures from RCSB when possible and
falls back to representative local SMOG-CHECK PDBs when downloads fail or when
official SMOG2 rejects a downloaded file.

Run from the repository root:

```bash
bash scripts/validate_real_pdb_panel.sh
bash scripts/validate_real_pdb_panel.sh --use-installed-smog3
bash scripts/validate_real_pdb_panel.sh --local-only --cases protein_ci2,dna
```

Generated inputs, outputs, and reports are written under ignored
`real_pdb_validation/`.
