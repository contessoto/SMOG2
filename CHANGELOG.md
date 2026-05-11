# Changelog

## 0.1.0a1

- Initial TestPyPI alpha release of `smog3`.
- Added Python-native SMOG3 runtime for the SMOG-CHECK-covered SMOG2 command
  surface.
- Verified original drop-in SMOG-CHECK harness tests 1 through 115 using a
  temporary SMOG3-backed `smog2` wrapper.
- Enforced zero-Perl SMOG3 normal runtime; Java remains supported for
  `SCM.jar` contact generation.
- Added two-stage parity workflows that compare original SMOG2 Docker baseline
  outputs with local SMOG3 candidate outputs.
- Added real-PDB panel validation helpers for comparing TestPyPI-installed or
  source-checkout SMOG3 against official SMOG2 Docker.
