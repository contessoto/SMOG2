from pathlib import Path

def _count_ndx_atoms(path: Path) -> int:
    return sum(
        len(line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("[")
    )


import smog3.cli as cli


CASES = [
    (53, "2ci2_v2.CB", ["-AAMATCH", "-contactMode", "shadow", "-contactParam", "10.0"]),
    (54, "1AKEapo_v2.CB", ["-AAMATCH", "-contactMode", "shadow", "-contactParam", "10.0"]),
    (55, "3IZH_v2.CB", ["-AAMATCH", "-contactMode", "shadow", "-contactParam", "10.0"]),
    (108, "1AKEapo_v3.BOND", ["-CABOND"]),
    (109, "glycans.BOND", ["-AABOND"]),
    (110, "RNA+protein", ["-AADIHE"]),
    (112, "RNA+protein", ["-AADIHE4"]),
]


def test_match_bond_dihe_group_runs_native_without_perl(monkeypatch, tmp_path: Path):
    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked: {argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)
    pdb_root = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files"

    for cid, stem, flags in CASES:
        pdb = pdb_root / f"{stem}.pdb"
        base = tmp_path / f"case_{cid}"
        top = base.with_suffix('.top')
        gro = base.with_suffix('.gro')
        ndx = base.with_suffix('.ndx')
        contacts = base.with_suffix('.contacts')

        monkeypatch.setattr(cli.sys, "argv", ["smog2", "-i", str(pdb), *flags, "-o", str(top), "-g", str(gro), "-n", str(ndx), "-s", str(contacts)])
        assert cli.smog2_main() == 0

        natoms = sum(1 for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM")))
        top_txt = top.read_text()
        assert "[ atoms ]" in top_txt and "[ bonds ]" in top_txt and "[ molecules ]" in top_txt
        assert int(gro.read_text().splitlines()[1].strip()) == natoms
        assert _count_ndx_atoms(ndx) == natoms

        crows = [ln for ln in contacts.read_text().splitlines() if ln.strip()]
        if "-contactMode" in flags:
            assert len(crows) > 0

    assert called["perl"] is False
