from pathlib import Path

def _count_ndx_atoms(path: Path) -> int:
    return sum(
        len(line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("[")
    )


import smog3.cli as cli
from smog3 import smog2_native


def _count_section(top_text: str, section: str) -> int:
    n = 0
    active = False
    for ln in top_text.splitlines():
        s = ln.strip()
        if s.startswith("[") and s.endswith("]"):
            active = s == f"[ {section} ]"
            continue
        if active and s and not s.startswith(";"):
            n += 1
    return n


def _read_contacts(path: Path):
    out = []
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith((";", "#")):
            continue
        out.append(tuple(s.split()))
    return out


def test_aa_test_template_uses_constant_fes_bond_length():
    template = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "templates" / "SBM_AA" / "AA-test.bif"
    atoms = [
        (1, "S1", "FES", 1, 0.0, 0.0, 0.0, "A"),
        (2, "FE1", "FES", 1, 2.2, 0.0, 0.0, "A"),
    ]
    attrs = smog2_native._template_atom_attributes(template)
    rules = smog2_native._template_bond_length_rules(template)

    assert smog2_native._distance_nm(atoms, 1, 2) == 0.22000000000000003
    assert smog2_native._template_bond_length_override(atoms, ["FES", "FES"], attrs, rules, 1, 2) == 0.21


def test_user_contact_and_2cg_cases_use_native_without_perl(monkeypatch, tmp_path: Path):
    # Case 50: 2ci2_v2 AA default-userC
    # Case 51: 2ci2_v2 AA default-gaussian-userC
    # Case 52: protein-RNA AA-2cg default
    pdb_root = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files"
    userc = pdb_root / "2ci2_v2.contacts"

    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked: {argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)

    cases = [
        ("2ci2_v2", ["-AA", "-userContacts", str(userc)], "AA"),
        ("2ci2_v2", ["-AAgaussian", "-userContacts", str(userc)], "AA"),
        ("protein-RNA", ["-AA2cg"], "AA2CG"),
    ]

    for stem, mode_args, atomtype in cases:
        pdb = pdb_root / f"{stem}.pdb"
        base = tmp_path / (stem.replace("+", "_").replace("-", "_") + "_" + atomtype)

        top = base.with_suffix(".top")
        gro = base.with_suffix(".gro")
        ndx = base.with_suffix(".ndx")
        contacts = base.with_suffix(".contacts")

        monkeypatch.setattr(
            cli.sys,
            "argv",
            ["smog2", "-i", str(pdb), *mode_args, "-o", str(top), "-g", str(gro), "-n", str(ndx), "-s", str(contacts)],
        )
        rc = cli.smog2_main()
        assert rc == 0

        top_text = top.read_text()
        pdb_atoms = sum(1 for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM")))
        pdb_res = len({(ln[21:22], ln[22:26].strip(), ln[26:27], ln[17:20].strip()) for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM"))})

        assert _count_section(top_text, "atoms") == pdb_atoms
        if atomtype == "AA2CG":
            assert _count_section(top_text, "bonds") > max(0, pdb_atoms - 1)
        else:
            assert _count_section(top_text, "bonds") == max(0, pdb_atoms - 1)
        assert _count_section(top_text, "molecules") == 1
        assert pdb_res >= 1
        assert (" NB_1 " in top_text) if atomtype == "AA2CG" else (f" {atomtype} " in top_text)

        ndx_atoms = _count_ndx_atoms(ndx)
        assert ndx_atoms == pdb_atoms
        assert int(gro.read_text().splitlines()[1].strip()) == pdb_atoms

        contact_rows = _read_contacts(contacts)
        if "-userContacts" in mode_args:
            src_rows = _read_contacts(userc)
            assert contact_rows == src_rows
            assert len(contact_rows) > 0
        elif atomtype == "AA2CG":
            assert len(contact_rows) > 0
        else:
            assert len(contact_rows) == 0

        for section in ("[ atoms ]", "[ bonds ]", "[ angles ]", "[ dihedrals ]", "[ pairs ]", "[ molecules ]"):
            assert section in top_text

    assert called["perl"] is False
