from pathlib import Path

from smog3.adjustpdb_native import main as adjust_main


def test_adjust_altloc_selection(tmp_path: Path):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text(
        """ATOM      1  CA AALA A   1      1.000   2.000   3.000  1.00  0.00           C
ATOM      2  CA BALA A   1      1.100   2.100   3.100  1.00  0.00           C
"""
    )
    rc = adjust_main(["-i", str(inp), "-o", str(out)])
    assert rc == 0
    txt = out.read_text()
    assert "BALA" not in txt


def test_adjust_remove_h_and_water(tmp_path: Path):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text(
        """ATOM      1  H1  ALA A   1      1.000   2.000   3.000  1.00  0.00           H
HETATM    2  O   HOH A   2      1.100   2.100   3.100  1.00  0.00           O
ATOM      3  CA  ALA A   3      1.200   2.200   3.200  1.00  0.00           C
"""
    )
    rc = adjust_main(["-i", str(inp), "-o", str(out), "-removeH", "-removewater"])
    assert rc == 0
    lines = [l for l in out.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert len(lines) == 1


def test_adjust_pdbnums_preserved(tmp_path: Path):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text(
        """ATOM     10  CA  ALA A  11      1.000   2.000   3.000  1.00  0.00           C
ATOM     20  CA  GLY A  21      1.100   2.100   3.100  1.00  0.00           C
"""
    )
    rc = adjust_main(["-i", str(inp), "-o", str(out), "-PDBnums"])
    assert rc == 0
    atom_lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
    assert atom_lines[0][6:11].strip() == "10"
    assert atom_lines[1][6:11].strip() == "20"
