# Changelog

## 0.1.0a3

- Prepared the next TestPyPI alpha after completing public tutorial
  model-generation parity.
- Verified `PYTHONPATH=src pytest -q tests`: 120 passed.
- Verified the original drop-in SMOG-CHECK harness: passed tests 1 through 115
  using the SMOG3-backed `smog2` wrapper.
- Verified public SMOG tutorial generation workflows: 22 PASS, 0 DIFF,
  0 SMOG2/SMOG3 errors, 0 manual-input model-generation cases, and
  2 simulation-only `NOT_GENERATION_TEST` cases.
- Verified an installed wheel outside the source repository can generate AA,
  CA, and OpenSMOG AA models for 2CI2 without invoking Perl, SMOG2, or Docker.
- Confirmed `SCM.jar` is bundled as package data under
  `smog3/data/tools/SCM.jar` for Java SCM contact generation.

## 0.1.0a2

- Prepared a new TestPyPI alpha build after the original `0.1.0a1` package was
  published and installed successfully.
- Added public SMOG tutorial/model-generation validation against official
  `smogserver/smog2:stable` Docker baselines, including strict comparison for
  representative AA, CA, OpenSMOG, ligand, nucleic-acid, disulfide, and
  user-contact cases.
- Improved `smog3 --help` so supported OpenSMOG flags are visible:
  `-OpenSMOG -OpenSMOGxml model.xml`.
- Enabled the real Java-SCM contact/topology path by default for normal AA
  model generation, fixing empty-contact output from direct installed
  `smog3 -AA` runs.
- Kept the original drop-in SMOG-CHECK harness validation green for tests 1
  through 115 using the SMOG3-backed `smog2` wrapper.
- Confirmed the SMOG3 runtime remains Python plus Java SCM only, with zero Perl
  in normal runtime.
- Updated README and validation docs with the installed TestPyPI workflow.

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
- Added public SMOG tutorial/model-generation validation helpers with a manifest
  for implemented, missing-input, and simulation-only tutorial categories.
