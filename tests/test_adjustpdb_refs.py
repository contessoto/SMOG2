from pathlib import Path
from io import StringIO

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


def test_adjust_map_applies_terminal_residue_names(tmp_path: Path):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    mapping = tmp_path / "map"
    mapping.write_text(
        """residue ASN C CA CB CG N ND2 O OD1
residue ASNN C CA CB CG N ND2 O OD1 %first
residue GLY C CA N O
residue GLYC C CA N O OXT %last
"""
    )
    inp.write_text(
        """ATOM      1  N   ASN A   1      1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ASN A   1      1.100   2.100   3.100  1.00  0.00           C
ATOM      3  N   GLY A   2      2.000   3.000   4.000  1.00  0.00           N
ATOM      4  OXT GLY A   2      2.100   3.100   4.100  1.00  0.00           O
"""
    )

    rc = adjust_main(["-i", str(inp), "-o", str(out), "-map", str(mapping)])

    assert rc == 0
    atom_lines = [line for line in out.read_text().splitlines() if line.startswith("ATOM")]
    assert atom_lines[0][17:21].strip() == "ASNN"
    assert atom_lines[-1][17:21].strip() == "GLYC"


def test_adjust_insertter_uses_prompt_answers(tmp_path: Path, monkeypatch):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text(
        """ATOM      1  CA  ALA A   1      1.000   2.000   3.000  1.00  0.00           C
ATOM      2  CA  ALA A   4      1.100   2.100   3.100  1.00  0.00           C
ATOM      3  CA  ALA A   9      1.200   2.200   3.200  1.00  0.00           C
"""
    )
    monkeypatch.setattr("sys.stdin", StringIO("n\ny\n"))

    rc = adjust_main(["-i", str(inp), "-o", str(out), "-insertTER"])

    assert rc == 0
    lines = out.read_text().splitlines()
    assert sum(1 for line in lines if line.startswith("TER")) == 1
    assert lines.index("TER") > 1


def test_adjust_preserves_existing_ter_and_resets_chain_numbering(tmp_path: Path):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text(
        """ATOM     10  CA  ALA A  11      1.000   2.000   3.000  1.00  0.00           C
TER
ATOM     20  CA  GLY B  21      1.100   2.100   3.100  1.00  0.00           C
"""
    )

    rc = adjust_main(["-i", str(inp), "-o", str(out)])

    assert rc == 0
    lines = out.read_text().splitlines()
    assert "TER" in lines
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    assert atom_lines[0][6:11].strip() == "1"
    assert atom_lines[0][22:26].strip() == "1"
    assert atom_lines[1][6:11].strip() == "1"
    assert atom_lines[1][22:26].strip() == "1"
