from pathlib import Path

from smog3.editgro_native import main as editgro_main


def _make_gro(path: Path):
    path.write_text(
        """h
2
    1RES     A1    1   0.000   0.000   0.000
    1RES     A2    2   1.000   1.000   1.000
1.00000 1.00000 1.00000
"""
    )


def test_editgro_center_and_pbc(tmp_path: Path):
    inp = tmp_path / "in.gro"
    out = tmp_path / "out.gro"
    _make_gro(inp)
    rc = editgro_main(["-g", str(inp), "-og", str(out), "-boxtype", "cubic", "-d", "1", "-c", "-pbc"])
    assert rc == 0
    lines = out.read_text().splitlines()
    assert lines[1].strip() == "2"
    box = [float(x) for x in lines[-1].split()]
    assert box[0] == box[1] == box[2]
