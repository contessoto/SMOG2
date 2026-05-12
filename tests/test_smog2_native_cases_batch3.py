from pathlib import Path

def _count_ndx_atoms(path: Path) -> int:
    return sum(
        len(line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("[")
    )


import smog3.cli as cli


def _pdb_counts(pdb_path: Path):
    atom_count = 0
    residues = set()
    for ln in pdb_path.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            atom_count += 1
            residues.add((ln[21:22], ln[22:26].strip(), ln[26:27], ln[17:20].strip()))
    return atom_count, len(residues)


def _top_section_count(top_text: str, section: str) -> int:
    lines = top_text.splitlines()
    in_section = False
    n = 0
    marker = f"[ {section} ]"
    for ln in lines:
        s = ln.strip()
        if s.startswith("[") and s.endswith("]"):
            in_section = s == marker
            continue
        if in_section and s and not s.startswith(";"):
            n += 1
    return n


def test_next_ten_direct_cases_use_native_without_perl(monkeypatch, tmp_path: Path):
    # smog2.testlist cases 10-19 (AA default)
    stems = [
        "1reschains_v2",
        "1AKEapo_v2",
        "2ci2_v2",
        "4gvy",
        "3IZH_v2",
        "terminaltest.BOND",
        "1F4N_v2",
        "3PTA",
        "tRNA",
        "tRNA.chop",
    ]

    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked for argv={argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)
    pdb_root = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files"

    for stem in stems:
        pdb = pdb_root / f"{stem}.pdb"
        dname = tmp_path / stem.replace(".", "_")

        monkeypatch.setattr(cli.sys, "argv", ["smog2", "-i", str(pdb), "-AA", "-dname", str(dname)])
        rc = cli.smog2_main()
        assert rc == 0

        top = Path(f"{dname}.top")
        gro = Path(f"{dname}.gro")
        ndx = Path(f"{dname}.ndx")
        contacts = Path(f"{dname}.contacts")

        top_text = top.read_text()
        gro_lines = gro.read_text().splitlines()
        ndx_atoms = _count_ndx_atoms(ndx)
        contact_lines = [ln for ln in contacts.read_text().splitlines() if ln.strip()]

        pdb_atoms, pdb_residues = _pdb_counts(pdb)
        top_atoms = _top_section_count(top_text, "atoms")
        top_bonds = _top_section_count(top_text, "bonds")
        top_molecules = _top_section_count(top_text, "molecules")

        assert int(gro_lines[1].strip()) == pdb_atoms
        assert top_atoms == pdb_atoms
        assert ndx_atoms == pdb_atoms
        assert top_bonds >= max(0, pdb_atoms - 1)
        assert top_molecules == 1
        assert pdb_residues >= 1
        assert len(contact_lines) > 0

        for section in (
            "[ defaults ]",
            "[ atomtypes ]",
            "[ moleculetype ]",
            "[ atoms ]",
            "[ bonds ]",
            "[ angles ]",
            "[ dihedrals ]",
            "[ pairs ]",
            "[ exclusions ]",
            "[ system ]",
            "[ molecules ]",
        ):
            assert section in top_text

    assert called["perl"] is False
