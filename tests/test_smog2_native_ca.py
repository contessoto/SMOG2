from pathlib import Path

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
