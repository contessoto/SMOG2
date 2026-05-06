from pathlib import Path

from smog3.smog2_native import main as smog2_main


def test_smog2_native_slice_on_first_testlist_case(tmp_path: Path):
    # first direct case in SMOG-CHECK testlist uses 1A01-AMP with AA default
    pdb = Path("SMOG-CHECK/share/PDB.files/1A01-AMP.pdb")
    cwd = Path.cwd()
    try:
        import os
        os.chdir(tmp_path)
        rc = smog2_main(["-i", str((cwd / pdb).resolve()), "-AA", "-dname", "slice"])
        assert rc == 0
        assert Path("slice.top").exists()
        assert Path("slice.gro").exists()
        assert Path("slice.ndx").exists()
        assert "[ atoms ]" in Path("slice.top").read_text()
    finally:
        os.chdir(cwd)
