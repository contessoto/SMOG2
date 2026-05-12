# Tutorial Command Audit

- Latest run: `validation/tutorials/runs/20260512-130713`
- Cases audited: `24`
- SMOG3 Perl invocations: `0`

| Metric | Count |
| --- | ---: |
| cases_with_executed_smog2_and_smog3_and_comparisons | 18 |
| cases_with_exact_public_tutorial_workflow_commands_executed | 8 |
| cases_executing_smog_adjustPDB | 8 |
| cases_with_tutorial_steps_using_smog_adjustPDB | 18 |
| cases_executing_removewater | 8 |
| cases_with_tutorial_steps_using_removewater | 15 |
| cases_executing_smog_ions | 2 |
| smog_ions_commands_executed | 4 |
| cases_with_tutorial_steps_using_smog_ions | 2 |
| cases_executing_custom_templates_maps_or_contacts | 8 |
| cases_with_tutorial_steps_using_custom_templates_maps_or_contacts | 14 |

| Status | Count |
| --- | ---: |
| PASS | 16 |
| DIFF | 1 |
| SMOG2_ERROR | 1 |
| SMOG3_ERROR | 0 |
| MISSING_INPUT | 0 |
| MISSING_DOWNLOAD | 0 |
| MANUAL_INPUT_REQUIRED | 4 |
| UNSUPPORTED_BY_SMOG3 | 0 |
| NOT_GENERATION_TEST | 2 |

## Cases

### standard_aa_ci2

- Tutorial: Standard all-atom model
- Source: https://smog-server.org/tutorials/OpenSMOG.AA/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA/steps.OpenSMOG.AA.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -AA -OpenSMOG -dname 2ci2.OpenSMOG.AA
python run.ci2.OpenSMOG.AA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/standard_aa_ci2/2ci2_v2.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_aa_ci2/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_aa_ci2/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_aa_ci2/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_aa_ci2/model.contacts -AA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/standard_aa_ci2/2ci2_v2.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/standard_aa_ci2/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/standard_aa_ci2/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/standard_aa_ci2/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/standard_aa_ci2/model.contacts -AA
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### standard_ca_ci2

- Tutorial: Standard C-alpha model
- Source: https://smog-server.org/tutorials/OpenSMOG.CA/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.CA/steps.OpenSMOG.CA.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -CA -OpenSMOG -dname 2ci2.OpenSMOG.CA
python run.ci2.OpenSMOG.CA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/standard_ca_ci2/2ci2_v2.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_ca_ci2/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_ca_ci2/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_ca_ci2/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/standard_ca_ci2/model.contacts -CA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/standard_ca_ci2/2ci2_v2.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/standard_ca_ci2/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/standard_ca_ci2/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/standard_ca_ci2/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/standard_ca_ci2/model.contacts -CA
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### opensmog_aa_ci2

- Tutorial: Standard all-atom model for OpenMM/OpenSMOG
- Source: https://smog-server.org/tutorials/OpenSMOG.AA/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA/steps.OpenSMOG.AA.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -AA -OpenSMOG -dname 2ci2.OpenSMOG.AA
python run.ci2.OpenSMOG.AA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/opensmog_aa_ci2/2ci2_v2.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_aa_ci2/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_aa_ci2/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_aa_ci2/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_aa_ci2/model.contacts -AA -OpenSMOG -OpenSMOGxml /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_aa_ci2/model.xml'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/opensmog_aa_ci2/2ci2_v2.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_aa_ci2/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_aa_ci2/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_aa_ci2/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_aa_ci2/model.contacts -AA -OpenSMOG -OpenSMOGxml validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_aa_ci2/model.xml
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`
- `model.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)

### opensmog_ca_ci2

- Tutorial: Standard C-alpha model for OpenMM/OpenSMOG
- Source: https://smog-server.org/tutorials/OpenSMOG.CA/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.CA/steps.OpenSMOG.CA.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -CA -OpenSMOG -dname 2ci2.OpenSMOG.CA
python run.ci2.OpenSMOG.CA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/opensmog_ca_ci2/2ci2_v2.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_ca_ci2/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_ca_ci2/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_ca_ci2/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_ca_ci2/model.contacts -CA -OpenSMOG -OpenSMOGxml /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/opensmog_ca_ci2/model.xml'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/opensmog_ca_ci2/2ci2_v2.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_ca_ci2/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_ca_ci2/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_ca_ci2/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_ca_ci2/model.contacts -CA -OpenSMOG -OpenSMOGxml validation/tutorials/runs/20260512-130713/smog3_candidate/opensmog_ca_ci2/model.xml
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`
- `model.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)

### aa_1a01_amp_ligand

- Tutorial: All-atom model with a novel ligand
- Source: https://smog-server.org/tutorials/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: none
- smog_adjustPDB executed / mentioned: `False` / `False`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `False` / `False`
- SMOG3 invoked Perl: `False`

SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/aa_1a01_amp_ligand/1A01-AMP.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_1a01_amp_ligand/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_1a01_amp_ligand/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_1a01_amp_ligand/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_1a01_amp_ligand/model.contacts -AA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/aa_1a01_amp_ligand/1A01-AMP.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/aa_1a01_amp_ligand/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/aa_1a01_amp_ligand/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/aa_1a01_amp_ligand/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/aa_1a01_amp_ligand/model.contacts -AA
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### aa_dna_terminal

- Tutorial: All-atom model for DNA terminal residues
- Source: https://smog-server.org/tutorials/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: none
- smog_adjustPDB executed / mentioned: `False` / `False`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `False` / `False`
- SMOG3 invoked Perl: `False`

SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/aa_dna_terminal/DNA.terminal.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_dna_terminal/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_dna_terminal/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_dna_terminal/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_dna_terminal/model.contacts -AA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/aa_dna_terminal/DNA.terminal.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/aa_dna_terminal/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/aa_dna_terminal/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/aa_dna_terminal/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/aa_dna_terminal/model.contacts -AA
```
Files compared:
- `model.top`: `PASS` (topology header metadata, tiny floating-point print ULPs, and dihedral +/-180 endpoint print convention)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### aa_trna

- Tutorial: All-atom RNA model
- Source: https://smog-server.org/tutorials/
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: none
- smog_adjustPDB executed / mentioned: `False` / `False`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `False` / `False`
- SMOG3 invoked Perl: `False`

SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/aa_trna/tRNA.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_trna/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_trna/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_trna/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_trna/model.contacts -AA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/aa_trna/tRNA.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/aa_trna/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/aa_trna/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/aa_trna/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/aa_trna/model.contacts -AA
```
Files compared:
- `model.top`: `PASS` (topology header metadata, tiny floating-point print ULPs, and dihedral +/-180 endpoint print convention)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### aa_disulfide_local

- Tutorial: All-atom model with a disulfide bond
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+disulfide
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+disulfide/steps.OpenSMOG.AA+disulfide.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `True`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2hi7.mod.pdb -o 2hi7.adjusted.pdb -insertTER
smog2 -i 2hi7.adjusted.pdb -AA -OpenSMOG -dname 2hi7.OpenSMOG.AA+disulfide -warn 2
python run.2hi7.OpenSMOG.AA+disulfide.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/aa_disulfide_local/terminaltest.BOND.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_disulfide_local/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_disulfide_local/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_disulfide_local/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_disulfide_local/model.contacts -AA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/aa_disulfide_local/terminaltest.BOND.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/aa_disulfide_local/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/aa_disulfide_local/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/aa_disulfide_local/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/aa_disulfide_local/model.contacts -AA
```
Files compared:
- `model.top`: `PASS` (topology header metadata, tiny floating-point print ULPs, and dihedral +/-180 endpoint print convention)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### ca_disulfide_local

- Tutorial: C-alpha model with a disulfide bond
- Source: https://smog-server.org/tutorials/OpenSMOG.CA+disulfide
- Status: `PASS`
- Downloaded assets used: 1
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.CA+disulfide/steps.OpenSMOG.CA+disulfide.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `True`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2hi7.mod.pdb -o 2hi7.adjusted.pdb -insertTER
smog2 -i 2hi7.adjusted.pdb -CA -OpenSMOG -dname 2hi7.OpenSMOG.CA+disulfide -warn 1
python run.2hi7.OpenSMOG.AA+disulfide.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/ca_disulfide_local/1AKEapo_v3.BOND.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/ca_disulfide_local/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/ca_disulfide_local/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/ca_disulfide_local/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/ca_disulfide_local/model.contacts -CA'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/ca_disulfide_local/1AKEapo_v3.BOND.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/ca_disulfide_local/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/ca_disulfide_local/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/ca_disulfide_local/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/ca_disulfide_local/model.contacts -CA
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS`

### aa_user_contacts_ci2

- Tutorial: All-atom model with multiple types of native contacts
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+multipleContactTypes
- Status: `PASS`
- Downloaded assets used: 2
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+multipleContactTypes/steps.OpenSMOG.AA+multipleContactTypes.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `True`
- OpenSMOG XML expected / compared: `False` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
smog2 -i 2ci2.adjusted.pdb -t SBM_AA-multipleContactTypes -opensmog -dname 2ci2.multipleContactTypes
python python run.multipleTypes.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir; smog2 -i /workdir/validation/tutorials/runs/20260512-130713/inputs/aa_user_contacts_ci2/2ci2_v2.pdb -SCMorig -keep4SCM -o /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_user_contacts_ci2/model.top -g /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_user_contacts_ci2/model.gro -n /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_user_contacts_ci2/model.ndx -s /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_user_contacts_ci2/model.contacts -AA -c /workdir/validation/tutorials/runs/20260512-130713/inputs/aa_user_contacts_ci2/2ci2_v2.contacts'
```
SMOG3 commands executed:
```text
python3 -m smog3.smog2_native -i validation/tutorials/runs/20260512-130713/inputs/aa_user_contacts_ci2/2ci2_v2.pdb -o validation/tutorials/runs/20260512-130713/smog3_candidate/aa_user_contacts_ci2/model.top -g validation/tutorials/runs/20260512-130713/smog3_candidate/aa_user_contacts_ci2/model.gro -n validation/tutorials/runs/20260512-130713/smog3_candidate/aa_user_contacts_ci2/model.ndx -s validation/tutorials/runs/20260512-130713/smog3_candidate/aa_user_contacts_ci2/model.contacts -AA -c validation/tutorials/runs/20260512-130713/inputs/aa_user_contacts_ci2/2ci2_v2.contacts
```
Files compared:
- `model.top`: `PASS` (topology header metadata before first section)
- `model.gro`: `PASS`
- `model.ndx`: `PASS`
- `model.contacts`: `PASS` (both files absent)

### aa_coulomb_electrostatics_tutorial

- Tutorial: All-atom model with Coulomb electrostatics
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+coulomb/
- Status: `PASS`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+coulomb/steps.OpenSMOG.AA+coulomb.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -t SBM_AA+coulomb -OpenSMOG -dname 2ci2.OpenSMOG.AA+coulomb
python run.ci2.OpenSMOG.AA+coulomb.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_coulomb_electrostatics_tutorial; smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_coulomb_electrostatics_tutorial; smog2 -i 2ci2.adjusted.pdb -t SBM_AA+coulomb -OpenSMOG -dname 2ci2.OpenSMOG.AA+coulomb'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
python3 -m smog3.smogcheck_dropin_smog2 -i 2ci2.adjusted.pdb -t SBM_AA+coulomb -OpenSMOG -dname 2ci2.OpenSMOG.AA+coulomb
```
Files compared:
- `2ci2.OpenSMOG.AA+coulomb.top`: `PASS` (topology header metadata before first section)
- `2ci2.OpenSMOG.AA+coulomb.gro`: `PASS`
- `2ci2.OpenSMOG.AA+coulomb.ndx`: `PASS`
- `2ci2.OpenSMOG.AA+coulomb.contacts`: `PASS`
- `2ci2.OpenSMOG.AA+coulomb.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)

### aa_debye_huckel_homogeneous_tutorial

- Tutorial: All-atom model with Debye-Huckel electrostatics and homogeneous nonbonded parameters
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+DebyeHuckel/
- Status: `PASS`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+DebyeHuckel/steps.OpenSMOG.AA+DebyeHuckel.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -t SBM_AA+DebyeHuckel -OpenSMOG -dname 2ci2.OpenSMOG.AA+DebyeHuckel
python run.ci2.OpenSMOG.AA+DebyeHuckel.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_debye_huckel_homogeneous_tutorial; smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_debye_huckel_homogeneous_tutorial; smog2 -i 2ci2.adjusted.pdb -t SBM_AA+DebyeHuckel -OpenSMOG -dname 2ci2.OpenSMOG.AA+DebyeHuckel'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
python3 -m smog3.smogcheck_dropin_smog2 -i 2ci2.adjusted.pdb -t SBM_AA+DebyeHuckel -OpenSMOG -dname 2ci2.OpenSMOG.AA+DebyeHuckel
```
Files compared:
- `2ci2.OpenSMOG.AA+DebyeHuckel.top`: `PASS` (topology header metadata before first section)
- `2ci2.OpenSMOG.AA+DebyeHuckel.gro`: `PASS`
- `2ci2.OpenSMOG.AA+DebyeHuckel.ndx`: `PASS`
- `2ci2.OpenSMOG.AA+DebyeHuckel.contacts`: `PASS`
- `2ci2.OpenSMOG.AA+DebyeHuckel.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)

### aa_explicit_ions_coulomb_tutorial

- Tutorial: All-atom model with explicit ions and Coulomb electrostatics
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+coulomb+ions/
- Status: `PASS`
- Downloaded assets used: 7
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+coulomb+ions/steps.OpenSMOG.AA+coulomb.ions.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `True` / `True`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
smog2 -i 2ci2.adjusted.pdb -t SBM_AA+coulomb+ions -OpenSMOG -dname 2ci2.OpenSMOG.AA+coulomb+ions -center -boxbuffer 5
smog_ions -t SBM_AA+coulomb+ions -f 2ci2.OpenSMOG.AA+coulomb+ions.top -g 2ci2.OpenSMOG.AA+coulomb+ions.gro -OpenSMOG 2ci2.OpenSMOG.AA+coulomb+ions.xml -of 2ci2.OpenSMOG.AA+coulomb+ions.MG.top -og 2ci2.OpenSMOG.AA+coulomb+ions.MG.gro -OpenSMOGout 2ci2.OpenSMOG.AA+coulomb+ions.MG.xml -ionnm MG
smog_ions -t SBM_AA+coulomb+ions -f 2ci2.OpenSMOG.AA+coulomb+ions.MG.top -g 2ci2.OpenSMOG.AA+coulomb+ions.MG.gro -OpenSMOG 2ci2.OpenSMOG.AA+coulomb+ions.MG.xml -of 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.top -og 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.gro -OpenSMOGout 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.xml -ionnm CL -ionn 10
python run.ci2.OpenSMOG.AA+coulomb+ions.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_coulomb_tutorial; smog_adjustPDB -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_coulomb_tutorial; smog2 -i 2ci2.adjusted.pdb -t SBM_AA+coulomb+ions -OpenSMOG -dname 2ci2.OpenSMOG.AA+coulomb+ions -center -boxbuffer 5'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_coulomb_tutorial; smog_ions -t SBM_AA+coulomb+ions -f 2ci2.OpenSMOG.AA+coulomb+ions.top -g 2ci2.OpenSMOG.AA+coulomb+ions.gro -OpenSMOG 2ci2.OpenSMOG.AA+coulomb+ions.xml -of 2ci2.OpenSMOG.AA+coulomb+ions.MG.top -og 2ci2.OpenSMOG.AA+coulomb+ions.MG.gro -OpenSMOGout 2ci2.OpenSMOG.AA+coulomb+ions.MG.xml -ionnm MG -ionn 2'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_coulomb_tutorial; smog_ions -t SBM_AA+coulomb+ions -f 2ci2.OpenSMOG.AA+coulomb+ions.MG.top -g 2ci2.OpenSMOG.AA+coulomb+ions.MG.gro -OpenSMOG 2ci2.OpenSMOG.AA+coulomb+ions.MG.xml -of 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.top -og 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.gro -OpenSMOGout 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.xml -ionnm CL -ionn 10'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -o 2ci2.adjusted.pdb -removewater
python3 -m smog3.smogcheck_dropin_smog2 -i 2ci2.adjusted.pdb -t SBM_AA+coulomb+ions -OpenSMOG -dname 2ci2.OpenSMOG.AA+coulomb+ions -center -boxbuffer 5
python3 -m smog3.ions_native -t SBM_AA+coulomb+ions -f 2ci2.OpenSMOG.AA+coulomb+ions.top -g 2ci2.OpenSMOG.AA+coulomb+ions.gro -OpenSMOG 2ci2.OpenSMOG.AA+coulomb+ions.xml -of 2ci2.OpenSMOG.AA+coulomb+ions.MG.top -og 2ci2.OpenSMOG.AA+coulomb+ions.MG.gro -OpenSMOGout 2ci2.OpenSMOG.AA+coulomb+ions.MG.xml -ionnm MG -ionn 2
python3 -m smog3.ions_native -t SBM_AA+coulomb+ions -f 2ci2.OpenSMOG.AA+coulomb+ions.MG.top -g 2ci2.OpenSMOG.AA+coulomb+ions.MG.gro -OpenSMOG 2ci2.OpenSMOG.AA+coulomb+ions.MG.xml -of 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.top -og 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.gro -OpenSMOGout 2ci2.OpenSMOG.AA+coulomb+ions.MGCL.xml -ionnm CL -ionn 10
```
Files compared:
- `2ci2.OpenSMOG.AA+coulomb+ions.top`: `PASS` (topology header metadata before first section)
- `2ci2.OpenSMOG.AA+coulomb+ions.gro`: `PASS` (smog_ions places ions stochastically; compared non-ion GRO records exactly, ion species/counts exactly, and candidate ion coordinates inside the same box)
- `2ci2.OpenSMOG.AA+coulomb+ions.ndx`: `PASS`
- `2ci2.OpenSMOG.AA+coulomb+ions.contacts`: `PASS`
- `2ci2.OpenSMOG.AA+coulomb+ions.xml`: `PASS` (OpenSMOG XML comments/whitespace and nonbond_param order)
- `2ci2.OpenSMOG.AA+coulomb+ions.MG.top`: `PASS` (topology header metadata and atomtypes whitespace rewritten by smog_ions)
- `2ci2.OpenSMOG.AA+coulomb+ions.MG.gro`: `PASS` (smog_ions places ions stochastically; compared non-ion GRO records exactly, ion species/counts exactly, and candidate ion coordinates inside the same box)
- `2ci2.OpenSMOG.AA+coulomb+ions.MG.xml`: `PASS` (OpenSMOG XML comments/whitespace and nonbond_param order)
- `2ci2.OpenSMOG.AA+coulomb+ions.MGCL.top`: `PASS` (topology header metadata and atomtypes whitespace rewritten by smog_ions)
- `2ci2.OpenSMOG.AA+coulomb+ions.MGCL.gro`: `PASS` (smog_ions places ions stochastically; compared non-ion GRO records exactly, ion species/counts exactly, and candidate ion coordinates inside the same box)
- `2ci2.OpenSMOG.AA+coulomb+ions.MGCL.xml`: `PASS` (OpenSMOG XML comments/whitespace and nonbond_param order)

### aa_explicit_ions_custom_potentials_tutorial

- Tutorial: All-atom model with explicit ions and custom effective potentials
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+custom+ions/
- Status: `DIFF`
- Downloaded assets used: 6
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+custom+ions/steps.OpenSMOG.AA+custom+ions.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `True` / `True`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -map AA_ions_Wang22.v1/AA_ions_Wang22.v1.map -removewater -o 2ci2.adjusted.pdb
smog2 -i 2ci2.adjusted.pdb -t AA_ions_Wang22.v1 -OpenSMOG -dname 2ci2.OpenSMOG.AA+custom+ions -center -boxbuffer 10
smog_ions -t AA_ions_Wang22.v1 -f 2ci2.OpenSMOG.AA+custom+ions.top -g 2ci2.OpenSMOG.AA+custom+ions.gro -OpenSMOG 2ci2.OpenSMOG.AA+custom+ions.xml -of 2ci2.OpenSMOG.AA+custom+ions.MG.top -og 2ci2.OpenSMOG.AA+custom+ions.MG.gro -OpenSMOGout 2ci2.OpenSMOG.AA+custom+ions.MG.xml -ionnm MG -ionn 2
smog_ions -t AA_ions_Wang22.v1 -f 2ci2.OpenSMOG.AA+custom+ions.MG.top -g 2ci2.OpenSMOG.AA+custom+ions.MG.gro -OpenSMOG 2ci2.OpenSMOG.AA+custom+ions.MG.xml -of 2ci2.OpenSMOG.AA+custom+ions.MGCL.top -og 2ci2.OpenSMOG.AA+custom+ions.MGCL.gro -OpenSMOGout 2ci2.OpenSMOG.AA+custom+ions.MGCL.xml -ionnm CL -ionn 10
python run.ci2.OpenSMOG.AA+custom+ions.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_custom_potentials_tutorial; smog_adjustPDB -i 2ci2.pdb -map AA_ions_Wang22.v1/AA_ions_Wang22.v1.map -removewater -o 2ci2.adjusted.pdb'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_custom_potentials_tutorial; smog2 -i 2ci2.adjusted.pdb -t AA_ions_Wang22.v1 -OpenSMOG -dname 2ci2.OpenSMOG.AA+custom+ions -center -boxbuffer 10'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_custom_potentials_tutorial; smog_ions -t AA_ions_Wang22.v1 -f 2ci2.OpenSMOG.AA+custom+ions.top -g 2ci2.OpenSMOG.AA+custom+ions.gro -OpenSMOG 2ci2.OpenSMOG.AA+custom+ions.xml -of 2ci2.OpenSMOG.AA+custom+ions.MG.top -og 2ci2.OpenSMOG.AA+custom+ions.MG.gro -OpenSMOGout 2ci2.OpenSMOG.AA+custom+ions.MG.xml -ionnm MG -ionn 2'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_explicit_ions_custom_potentials_tutorial; smog_ions -t AA_ions_Wang22.v1 -f 2ci2.OpenSMOG.AA+custom+ions.MG.top -g 2ci2.OpenSMOG.AA+custom+ions.MG.gro -OpenSMOG 2ci2.OpenSMOG.AA+custom+ions.MG.xml -of 2ci2.OpenSMOG.AA+custom+ions.MGCL.top -og 2ci2.OpenSMOG.AA+custom+ions.MGCL.gro -OpenSMOGout 2ci2.OpenSMOG.AA+custom+ions.MGCL.xml -ionnm CL -ionn 10'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -map AA_ions_Wang22.v1/AA_ions_Wang22.v1.map -removewater -o 2ci2.adjusted.pdb
python3 -m smog3.smogcheck_dropin_smog2 -i 2ci2.adjusted.pdb -t AA_ions_Wang22.v1 -OpenSMOG -dname 2ci2.OpenSMOG.AA+custom+ions -center -boxbuffer 10
python3 -m smog3.ions_native -t AA_ions_Wang22.v1 -f 2ci2.OpenSMOG.AA+custom+ions.top -g 2ci2.OpenSMOG.AA+custom+ions.gro -OpenSMOG 2ci2.OpenSMOG.AA+custom+ions.xml -of 2ci2.OpenSMOG.AA+custom+ions.MG.top -og 2ci2.OpenSMOG.AA+custom+ions.MG.gro -OpenSMOGout 2ci2.OpenSMOG.AA+custom+ions.MG.xml -ionnm MG -ionn 2
python3 -m smog3.ions_native -t AA_ions_Wang22.v1 -f 2ci2.OpenSMOG.AA+custom+ions.MG.top -g 2ci2.OpenSMOG.AA+custom+ions.MG.gro -OpenSMOG 2ci2.OpenSMOG.AA+custom+ions.MG.xml -of 2ci2.OpenSMOG.AA+custom+ions.MGCL.top -og 2ci2.OpenSMOG.AA+custom+ions.MGCL.gro -OpenSMOGout 2ci2.OpenSMOG.AA+custom+ions.MGCL.xml -ionnm CL -ionn 10
```
Files compared:
- `2ci2.OpenSMOG.AA+custom+ions.top`: `DIFF`
- `2ci2.OpenSMOG.AA+custom+ions.gro`: `DIFF` (non-ion GRO atom records differ)
- `2ci2.OpenSMOG.AA+custom+ions.ndx`: `PASS`
- `2ci2.OpenSMOG.AA+custom+ions.contacts`: `PASS`
- `2ci2.OpenSMOG.AA+custom+ions.xml`: `DIFF` (OpenSMOG XML semantic content differs)
- `2ci2.OpenSMOG.AA+custom+ions.MG.top`: `DIFF`
- `2ci2.OpenSMOG.AA+custom+ions.MG.gro`: `DIFF` (non-ion GRO atom records differ)
- `2ci2.OpenSMOG.AA+custom+ions.MG.xml`: `DIFF` (OpenSMOG XML semantic content differs)
- `2ci2.OpenSMOG.AA+custom+ions.MGCL.top`: `DIFF`
- `2ci2.OpenSMOG.AA+custom+ions.MGCL.gro`: `DIFF` (non-ion GRO atom records differ)
- `2ci2.OpenSMOG.AA+custom+ions.MGCL.xml`: `DIFF` (OpenSMOG XML semantic content differs)
Remaining reason/blocker: nual.
                                            
           For questions regarding this script, contact info@smog-server.org
****************************************************************************************
Custom nonbond function found in templates.
Will add ion nonbond terms to XML file, if needed.

Will read definitions of ion parameters in AA_ions_Wang22.v1/AA_ions_Wang22.v1.ions.def

Will add ions to 2ci2.OpenSMOG.AA+custom+ions.MG.top and 2ci2.OpenSMOG.AA+custom+ions.MG.gro
Will try to add 10 CL ions with charge -1 and mass 1
Will use C12=5.96046e-09 E*nm^12 and C6=0 E*nm^6
Will place ions at least 0.5 nm from the existing system
Will write the new force field and configuration to 2ci2.OpenSMOG.AA+custom+ions.MGCL.top, 2ci2.OpenSMOG.AA+custom+ions.MGCL.xml and 2ci2.OpenSMOG.AA+custom+ions.MGCL.gro

Rectangular box detected
box vectors:
	22.9164	0	0
	0	23.4903	0
	0	0	23.0271

preparing to add ions to gro file...
adding ions to gro file...
adding ions to top file...
atomt

### aa_modified_residues_tutorial

- Tutorial: All-atom model with non-canonical/modified residues
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+newResidue
- Status: `MANUAL_INPUT_REQUIRED`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+newResidue/steps.OpenSMOG.AA+newResidue.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `True`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 3eao.noligands.pdb -map sbmMapExact+SEC -removewater -o 3eao.noligands.adjusted.pdb
smog2 -i 3eao.noligands.adjusted.pdb -t SBM_AA+SEC -dname 3eao.noligands
python run.3eao.OpenSMOG.AA.py
```
Remaining reason/blocker: Public tutorial assets are fetched, but this tutorial uses custom template/script steps that need explicit SMOG3 CLI mapping before automated output comparison is trustworthy.

### aa_glycans_tutorial

- Tutorial: All-atom model with glycans
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+glycans
- Status: `MANUAL_INPUT_REQUIRED`
- Downloaded assets used: 2
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+glycans/steps.OpenSMOG.AA+glycans.txt
- smog_adjustPDB executed / mentioned: `False` / `False`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `True`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog2 -i prefusion_with_glycans.pdb -t AA_glycans_Dodero21.v1 -OpenSMOG -dname spike.OpenSMOG.AA+glycans
python run.glycans.py
```
Remaining reason/blocker: Public tutorial assets are fetched, but this tutorial uses custom template/script steps that need explicit SMOG3 CLI mapping before automated output comparison is trustworthy.

### aa_custom_contact_potential_tutorial

- Tutorial: All-atom model with a custom contact potential
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+customContacts
- Status: `PASS`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+customContacts/steps.OpenSMOG.AA+customContacts.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
smog2 -t SBM_AA+customContacts -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.AA+customContacts
python run.ci2.OpenSMOG.AA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_custom_contact_potential_tutorial; smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_custom_contact_potential_tutorial; smog2 -t SBM_AA+customContacts -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.AA+customContacts'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
python3 -m smog3.smogcheck_dropin_smog2 -t SBM_AA+customContacts -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.AA+customContacts
```
Files compared:
- `2ci2.OpenSMOG.AA+customContacts.top`: `PASS` (topology header metadata before first section)
- `2ci2.OpenSMOG.AA+customContacts.gro`: `PASS`
- `2ci2.OpenSMOG.AA+customContacts.ndx`: `PASS`
- `2ci2.OpenSMOG.AA+customContacts.contacts`: `PASS`
- `2ci2.OpenSMOG.AA+customContacts.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)

### aa_custom_nonbonded_heterogeneous_tutorial

- Tutorial: All-atom model with a custom nonbonded potential with heterogeneous parameters
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+customNonbonded
- Status: `PASS`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+customNonbonded/steps.OpenSMOG.AA+customNonbonded.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
smog2 -t SBM_AA+customNonbonded -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.AA+customNonbonded
python run.ci2.OpenSMOG.AA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_custom_nonbonded_heterogeneous_tutorial; smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/aa_custom_nonbonded_heterogeneous_tutorial; smog2 -t SBM_AA+customNonbonded -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.AA+customNonbonded'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
python3 -m smog3.smogcheck_dropin_smog2 -t SBM_AA+customNonbonded -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.AA+customNonbonded
```
Files compared:
- `2ci2.OpenSMOG.AA+customNonbonded.top`: `PASS` (topology header metadata before first section)
- `2ci2.OpenSMOG.AA+customNonbonded.gro`: `PASS`
- `2ci2.OpenSMOG.AA+customNonbonded.ndx`: `PASS`
- `2ci2.OpenSMOG.AA+customNonbonded.contacts`: `PASS`
- `2ci2.OpenSMOG.AA+customNonbonded.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)

### ca_custom_nonbonded_tutorial

- Tutorial: C-alpha model with custom nonbonded interactions
- Source: https://smog-server.org/tutorials/OpenSMOG.CA+customNonbonded
- Status: `SMOG2_ERROR`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.CA+customNonbonded/steps.OpenSMOG.CA+customNonbonded.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
smog2 -t SBM_AA -tCG SBM_CA+customNonbonded -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.CA+customNonbonded
python run.ci2.OpenSMOG.AA.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/ca_custom_nonbonded_tutorial; smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/ca_custom_nonbonded_tutorial; smog2 -t SBM_AA -tCG SBM_CA+customNonbonded -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.CA+customNonbonded'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
python3 -m smog3.smogcheck_dropin_smog2 -t SBM_AA -tCG SBM_CA+customNonbonded -i 2ci2.adjusted.pdb -opensmog -dname 2ci2.OpenSMOG.CA+customNonbonded
```
Files compared:
- `2ci2.OpenSMOG.CA+customNonbonded.top`: `DIFF` (missing file)
- `2ci2.OpenSMOG.CA+customNonbonded.gro`: `DIFF`
- `2ci2.OpenSMOG.CA+customNonbonded.ndx`: `DIFF`
- `2ci2.OpenSMOG.CA+customNonbonded.contacts`: `DIFF`
- `2ci2.OpenSMOG.CA+customNonbonded.xml`: `DIFF` (missing file)
Remaining reason/blocker: bonded/CA+customContacts.extras
Found .nb file SBM_CA+customNonbonded/CA+customContacts.nb
Found .bif file SBM_CA+customNonbonded/CA+customContacts.bif
Found .sif file SBM_CA+customNonbonded/CA+customContacts.sif
Found .b file SBM_CA+customNonbonded/CA+customContacts.b
Validating XML format of SBM_CA+customNonbonded/CA+customContacts.bif
Validating XML format of SBM_CA+customNonbonded/CA+customContacts.sif


FATAL ERROR:  Failed at validating SBM_CA+customNonbonded/CA+customContacts.sif: Illegal attribute 'nbcutoff' found in <CustomNonBonded>.
 
This is due to an XML formatting issue.  The most common issue is that an element is missing a tagline. For example, something like this may appear in your file,
  <child>
    <subchild>.....
  </child>

whereas the following would be appropriate:
  <child>
    <subchild>.....

    </subchild>
  </child>

For more information about specific errors, you can check the FAQ page on smog-server.org,
the SMOG2 manual, or you can email us at info@smog

### large_system_fragment_tutorial

- Tutorial: Fragment of a large system using the all-atom model
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+largeFragment
- Status: `MANUAL_INPUT_REQUIRED`
- Downloaded assets used: 3
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+largeFragment/steps.OpenSMOG.AA.largeFragment.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 7k00.modified.pdb -insertTER -removewater -o 7k00.adjusted.pdb
smog2 -i 7k00.adjusted.pdb -AA -opensmog -dname 7k00.OpenSMOG.AA
smog_extract -f 7k00.OpenSMOG.AA.top -g 7k00.OpenSMOG.AA.gro -n keep.ndx -of ribosome.fragment.top -og ribosome.fragment.gro -OpenSMOG 7k00.OpenSMOG.AA.xml -OpenSMOGout ribosome.fragment.xml
python run.7k00.py
python run.fragment.py
```
Remaining reason/blocker: Public tutorial assets are fetched, but this tutorial uses custom template/script steps that need explicit SMOG3 CLI mapping before automated output comparison is trustworthy.

### checkpoint_many_segment_workflow

- Tutorial: Performing simulations with many segments via checkpoint files
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+checkpoints
- Status: `NOT_GENERATION_TEST`
- Downloaded assets used: 1
- Steps files: none
- smog_adjustPDB executed / mentioned: `False` / `False`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `False` / `False`
- SMOG3 invoked Perl: `False`


### steepest_descent_minimization

- Tutorial: Steepest descent minimization
- Source: https://smog-server.org/tutorials/OpenSMOG+GradientDescent
- Status: `NOT_GENERATION_TEST`
- Downloaded assets used: 1
- Steps files: none
- smog_adjustPDB executed / mentioned: `False` / `False`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `False`
- OpenSMOG XML expected / compared: `False` / `False`
- SMOG3 invoked Perl: `False`


### public_aa_novel_ligand_downloaded

- Tutorial: Public all-atom model with a novel ligand
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+novelLigand
- Status: `MANUAL_INPUT_REQUIRED`
- Downloaded assets used: 4
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+novelLigand/steps.OpenSMOG.AA+novelLigand.txt
- smog_adjustPDB executed / mentioned: `False` / `True`
- -removewater executed / mentioned: `False` / `False`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `False` / `True`
- OpenSMOG XML expected / compared: `True` / `False`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 3p63.pdb -o 3p63.adjusted.pdb -map sbmMapExact+FES -insertTER
smog2 -i 3p63.adjusted.pdb -t SBM_AA+novelLigand -OpenSMOG -dname 3p63.OpenSMOG.AA+novelLigand
python run.ci2.OpenSMOG.AA.py
```
Remaining reason/blocker: Public tutorial assets are fetched, but this tutorial uses custom template/script steps that need explicit SMOG3 CLI mapping before automated output comparison is trustworthy.

### public_aa_multiple_contact_types_downloaded

- Tutorial: Public all-atom model with multiple native-contact types
- Source: https://smog-server.org/tutorials/OpenSMOG.AA+multipleContactTypes
- Status: `PASS`
- Downloaded assets used: 3
- Steps files: validation/tutorials/assets/downloads/smog-server.org/tutorials/OpenSMOG.AA+multipleContactTypes/steps.OpenSMOG.AA+multipleContactTypes.txt
- smog_adjustPDB executed / mentioned: `True` / `True`
- -removewater executed / mentioned: `True` / `True`
- smog_ions executed / mentioned: `False` / `False`
- custom templates/maps/contact files executed / mentioned: `True` / `True`
- OpenSMOG XML expected / compared: `True` / `True`
- SMOG3 invoked Perl: `False`

Public tutorial commands found:
```text
smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
smog2 -i 2ci2.adjusted.pdb -t SBM_AA-multipleContactTypes -opensmog -dname 2ci2.multipleContactTypes
python python run.multipleTypes.py
```
SMOG2 Docker commands executed:
```text
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/public_aa_multiple_contact_types_downloaded; smog_adjustPDB -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb'
docker run --rm -v /Users/vc18/Work/Dev/SMOG2:/workdir smogserver/smog2:stable bash -lc 'set -euo pipefail; cd /workdir/validation/tutorials/runs/20260512-130713/smog2_baseline/public_aa_multiple_contact_types_downloaded; smog2 -i 2ci2.adjusted.pdb -t SBM_AA-multipleContactTypes -opensmog -dname 2ci2.multipleContactTypes'
```
SMOG3 commands executed:
```text
python3 -c 'from smog3.adjustpdb_native import main; import sys; raise SystemExit(main(sys.argv[1:]))' -i 2ci2.pdb -removewater -o 2ci2.adjusted.pdb
python3 -m smog3.smogcheck_dropin_smog2 -i 2ci2.adjusted.pdb -t SBM_AA-multipleContactTypes -opensmog -dname 2ci2.multipleContactTypes
```
Files compared:
- `2ci2.multipleContactTypes.top`: `PASS` (topology header metadata before first section)
- `2ci2.multipleContactTypes.gro`: `PASS`
- `2ci2.multipleContactTypes.ndx`: `PASS`
- `2ci2.multipleContactTypes.contacts`: `PASS`
- `2ci2.multipleContactTypes.xml`: `PASS` (OpenSMOG XML generated comment metadata before root element)
