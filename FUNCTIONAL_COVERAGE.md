# SMOG3 Functional Coverage Report

## Command coverage status

### smog2
- **Status:** partially native
- **Execution path:** `smog3.cli.smog2_main -> smog3.smog2_native.main (AA/CA default slice, including SMOG-CHECK direct AA cases 1-25, user-contact/AA-2cg cases 50-52, and OpenSMOG/g96 slice cases, plus shadow/cutoff parameterized contact-mode slice (cases 56-92)) with loud native errors for unsupported flag sets (optional legacy Perl fallback only via SMOG3_LEGACY_PERL_FALLBACK=1)`
- **Tests:** `tests/test_compat_wrappers.py`, `tests/test_smog2_native_slice.py`, `tests/test_smog2_native_cases.py`, `tests/test_smog2_native_cases_batch2.py`, `tests/test_smog2_native_cases_batch3.py`, `tests/test_smog2_native_cases_batch4.py`, `tests/test_smog2_native_cases_batch5.py`, `tests/test_smog2_native_feature_groups.py`, `tests/test_smog2_native_shadow_cutoff_group.py`, `tests/test_smog2_native_match_bond_dihe.py`, `tests/test_parity_direct.py`
- **Reference parity:** mixed (exact for selected legacy slices; structural/parity-algorithmic for newer native groups)
- **Remaining for full SMOG-CHECK:** close remaining partial/delegated template-map, freecoor/interactive, and advanced force-field parity gaps.

### smog_tablegen
- **Status:** fully native
- **Execution path:** `smog3.cli.smog_tablegen_main -> smog3.tablegen.main`
- **Tests:** `tests/test_tablegen_parity.py`, `tests/test_tablegen_refs.py`
- **Reference parity:** yes (`SMOG-CHECK/share/refs/table*.xvg`)
- **Remaining:** broader option-matrix parity.

### smog_scale-energies
- **Status:** partially native
- **Execution path:** `smog3.cli.smog_scale_energies_main -> smog3.scale_energies_native.main`
- **Tests:** `tests/test_scale_energies_native.py`
- **Reference parity:** synthetic parity only
- **Remaining:** full directive/comment/include parity and interactive group-selection parity.

### smog_extract
- **Status:** partially native
- **Execution path:** `smog3.cli.smog_extract_main -> smog3.extract_native.main`
- **Tests:** `tests/test_extract_native.py`
- **Reference parity:** synthetic parity only
- **Remaining:** full OpenSMOG semantics and full restraints/interaction-removal parity details.

### smog_adjustPDB
- **Status:** partially native
- **Execution path:** `smog3.cli.smog_adjustpdb_main -> smog3.adjustpdb_native.main`
- **Tests:** `tests/test_adjustpdb_refs.py` (logic-based synthetic cases)
- **Reference parity:** synthetic parity only
- **Remaining:** exact map/legacy semantics and byte-level parity across full SMOG-CHECK matrix.

### smog_editgro
- **Status:** partially native
- **Execution path:** `smog3.cli.smog_editgro_main -> smog3.editgro_native.main`
- **Tests:** `tests/test_editgro_native.py`
- **Reference parity:** synthetic parity only
- **Remaining:** strict parity vs `SMOG-CHECK/share/editgrorefs`.

### smog_ions
- **Status:** partially native
- **Execution path:** `smog3.cli.smog_ions_main -> smog3.ions_native.main`
- **Tests:** `tests/test_ions_native.py`
- **Reference parity:** synthetic parity only
- **Remaining:** full topology insertion ordering/format and OpenSMOG-aware behavior parity.

### smog_modifyXML
- **Status:** partially native
- **Execution path:** `smog3.cli.smog_modifyxml_main -> smog3.modifyxml_native.main`
- **Tests:** `tests/test_modifyxml_native.py`
- **Reference parity:** synthetic parity only
- **Remaining:** full interactive settings-driven behavior and XML schema-level parity with Perl implementation.


### smog3-parity-direct
- **Status:** available
- **Execution path:** `smog3.parity_direct.main` runs baseline Perl SMOG2 and Python-native SMOG3 in isolated temp dirs and compares `.top/.gro/.ndx/.contacts` (+`.xml` for OpenSMOG cases).
- **Default parity set:** cases `1,21,41,50,56,94`
- **Note:** if Perl deps are missing, baseline is skipped with explicit reason (non-zero exit).
