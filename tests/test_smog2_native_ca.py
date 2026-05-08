from pathlib import Path
from types import SimpleNamespace

import pytest

from smog3 import smog2_native
from smog3.smog2_native import main as smog2_main


def test_smog2_native_slice_supports_ca_model(tmp_path: Path):
    pdb = Path('SMOG-CHECK/share/PDB.files/2ci2_v2.pdb').resolve()
    rc = smog2_main(["-i", str(pdb), "-CA", "-dname", str(tmp_path / "ca_case")])
    assert rc == 0

    top = Path(f"{tmp_path / 'ca_case'}.top")
    assert top.exists()
    text = top.read_text()
    assert "[ atoms ]" in text
    assert " CA " in text


def test_ca_mapping_preserves_residue_numbering_chain_groups_and_terminal_name():
    pdb = Path('SMOG-CHECK/share/PDB.files/1reschains_v2.pdb').resolve()
    atoms = smog2_native._parse_pdb_atoms(pdb)
    ca_atoms = smog2_native._coarse_grain_ca_atoms(atoms)

    assert len(ca_atoms) == 210
    assert ca_atoms[0][:4] == (1, "CA", "MET", 1)
    assert ca_atoms[-1][:4] == (214, "CA", "GLYT", 214)
    assert [len(vals) for _chain, vals in smog2_native._chain_atom_groups(ca_atoms)] == [50, 1, 34, 1, 124]


def test_ca_bond_records_extend_topology_and_filter_bonded_contacts():
    pdb = Path('SMOG-CHECK/share/PDB.files/1AKEapo_v3.BOND.pdb').resolve()
    atoms = smog2_native._parse_pdb_atoms(pdb)
    ca_atoms = smog2_native._coarse_grain_ca_atoms(atoms)
    extra_bonds = smog2_native._ca_user_bonds_from_pdb(pdb, ca_atoms)

    assert extra_bonds == [(1, 2)]
    bonds, angles, dihedrals = smog2_native._ca_bonded_sections(ca_atoms, extra_bonds)
    assert (1, 2) in bonds
    assert (1, 2, 3) in angles
    assert any(row[:4] == (1, 2, 3, 4) for row in dihedrals)

    contacts = [
        (1, 1, ("2", "2")),
        (1, 1, ("2", "3")),
        (1, 1, ("2", "4")),
        (1, 1, ("2", "5")),
    ]
    assert smog2_native._filter_ca_contacts_by_bonded_exclusions(ca_atoms, contacts, extra_bonds) == [
        (1, 1, ("2", "5")),
    ]


def test_ca_scm_contacts_use_java_raw_atom_contacts_and_map_to_residues(monkeypatch, tmp_path: Path):
    source_atoms = [
        (1, "N", "ALA", 1, 0.0, 0.0, 0.0, "A:1"),
        (2, "CA", "ALA", 1, 1.0, 0.0, 0.0, "A:1"),
        (3, "N", "GLY", 5, 6.0, 0.0, 0.0, "A:1"),
        (4, "CA", "GLY", 5, 7.0, 0.0, 0.0, "A:1"),
    ]
    ca_atoms = smog2_native._coarse_grain_ca_atoms(source_atoms)
    calls = []

    def fake_run(argv, capture_output, text):
        calls.append(argv)
        assert Path(argv[0]).name != "perl"
        raw = Path(argv[argv.index("-o") + 1])
        raw.write_text("1 1 1 3\n1 2 1 4\n")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(smog2_native.shutil, "which", lambda name: f"/usr/bin/{name}" if name == "java" else None)
    monkeypatch.setattr(smog2_native.subprocess, "run", fake_run)

    contacts = smog2_native._generate_ca_contacts_with_scm(
        tmp_path / "model.gro",
        tmp_path / "model.top",
        tmp_path / "model.contacts",
        ca_atoms,
        source_atoms,
    )

    assert contacts == [(1, 1, ("1", "5"))]
    assert calls and calls[0][0] == "/usr/bin/java"
    assert calls[0][calls[0].index("-o") + 1].endswith("model.contacts.ShadowOutput")
    assert (tmp_path / "model.top4SCM.top").exists()
    gro4scm = (tmp_path / "model.gro4SCM.gro").read_text().splitlines()
    assert int(gro4scm[1]) == 2


def test_case41_ca_contacts_and_topology_when_scm_jar_available(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    if not (root / "src" / "tools" / "SCM.jar").exists():
        pytest.skip("SCM.jar is not available")
    if smog2_native.shutil.which("java") is None:
        pytest.skip("java is not available")

    pdb = root / "SMOG-CHECK" / "share" / "PDB.files" / "1reschains_v2.pdb"
    base = tmp_path / "case41"
    rc = smog2_main([
        "-i", str(pdb),
        "-CA",
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])

    assert rc == 0
    assert int(base.with_suffix(".gro").read_text().splitlines()[1]) == 210
    assert len([line for line in base.with_suffix(".contacts").read_text().splitlines() if line.strip()]) == 613
    assert base.with_suffix(".contacts").read_text().splitlines()[:3] == [
        "1 1 1 26",
        "1 1 3 79",
        "1 1 3 80",
    ]

    counts = {}
    current = None
    for raw in base.with_suffix(".top").read_text().splitlines():
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            counts[current] = 0
        elif current and s and not s.startswith(";"):
            counts[current] += 1
    assert counts["atoms"] == 210
    assert counts["bonds"] == 205
    assert counts["angles"] == 202
    assert counts["dihedrals"] == 398
    assert counts["pairs"] == 613
    assert counts["exclusions"] == 613
