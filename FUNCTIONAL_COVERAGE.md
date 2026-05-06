# SMOG3 Functional Coverage Report

## Command coverage status

### smog2
- **Status:** partially native
- **Execution path:** `smog3.cli.smog2_main -> smog3.smog2_native.main (AA slice) with fallback to perl for unsupported flag sets`
- **Tests:** `tests/test_compat_wrappers.py`, `tests/test_smog2_native_slice.py`
- **Reference parity:** not yet
- **Remaining for full SMOG-CHECK:** expand from AA default vertical slice to additional SMOG-CHECK cases and full template/contact semantics.

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
