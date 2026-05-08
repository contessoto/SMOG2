from pathlib import Path

def _count_ndx_atoms(path: Path) -> int:
    return sum(
        len(line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("[")
    )


import smog3.cli as cli


def _count_top_atoms(top_text: str) -> int:
    lines = top_text.splitlines()
    in_atoms = False
    n = 0
    for ln in lines:
        s = ln.strip()
        if s.startswith("[") and s.endswith("]"):
            in_atoms = s.strip("[] ") == "atoms"
            continue
        if in_atoms and s and not s.startswith(";"):
            n += 1
    return n


def test_next_five_direct_cases_use_native_without_perl(monkeypatch, tmp_path: Path):
    # smog2.testlist cases 5-9 (AA default)
    stems = [
        "2FP4-GDP",
        "2FP4-GTP",
        "2FP4-GDP.largeformat",
        "DNA.terminal.BMG",
        "DNA.terminal",
    ]

    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked for argv={argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)

    pdb_root = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files"

    for stem in stems:
        dname = tmp_path / stem.replace(".", "_")
        monkeypatch.setattr(cli.sys, "argv", ["smog2", "-i", str(pdb_root / f"{stem}.pdb"), "-AA", "-dname", str(dname)])
        rc = cli.smog2_main()
        assert rc == 0

        top = Path(f"{dname}.top")
        gro = Path(f"{dname}.gro")
        ndx = Path(f"{dname}.ndx")
        contacts = Path(f"{dname}.contacts")

        assert top.exists() and gro.exists() and ndx.exists() and contacts.exists()

        ttxt = top.read_text()
        natoms_top = _count_top_atoms(ttxt)
        natoms_gro = int(gro.read_text().splitlines()[1].strip())
        natoms_ndx = _count_ndx_atoms(ndx)

        for section in ("[ atoms ]", "[ bonds ]", "[ angles ]", "[ dihedrals ]", "[ pairs ]", "[ molecules ]"):
            assert section in ttxt

        assert natoms_top == natoms_gro == natoms_ndx
        assert contacts.read_text() == ""

    assert called["perl"] is False
