# Basic All-Atom Example

This example builds an all-atom SMOG3 model for 2CI2.

Run:

```bash
bash run.sh
```

The script downloads `2CI2.pdb` when possible and falls back to the local
SMOG-CHECK copy when run from a source checkout.

Generated files:

- `model.top`
- `model.gro`
- `model.ndx`
- `model.contacts`

Generated outputs are ignored and should not be committed.
