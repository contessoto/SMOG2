from pathlib import Path

def _count_ndx_atoms(path: Path) -> int:
    return sum(
        len(line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("[")
    )


from smog3.smog2_native import main as smog2_main


PDBS = [
    "1A01-ADP.pdb",
    "1A01-ATP.pdb",
    "1A01-ADP+RNA.pdb",
]


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


def test_next_three_direct_cases_run_natively(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    pdb_root = root / "SMOG-CHECK" / "share" / "PDB.files"

    for stem in PDBS:
        dname = stem.replace('.pdb', '').replace('+', '_')
        rc = smog2_main(["-i", str(pdb_root / stem), "-AA", "-dname", str(tmp_path / dname)])
        assert rc == 0

        top = Path(f"{tmp_path / dname}.top")
        gro = Path(f"{tmp_path / dname}.gro")
        ndx = Path(f"{tmp_path / dname}.ndx")
        contacts = Path(f"{tmp_path / dname}.contacts")

        ttxt = top.read_text()
        gtxt = gro.read_text().splitlines()
        ntx = ndx.read_text()
        natoms_gro = int(gtxt[1].strip())
        natoms_top = _count_top_atoms(ttxt)
        natoms_ndx = _count_ndx_atoms(ndx)

        assert "[ defaults ]" in ttxt
        assert "[ atomtypes ]" in ttxt
        assert "[ moleculetype ]" in ttxt
        assert "[ atoms ]" in ttxt
        assert "[ bonds ]" in ttxt
        assert "[ molecules ]" in ttxt
        assert natoms_gro == natoms_top == natoms_ndx
        assert contacts.read_text() == ""
