from pathlib import Path

from smog3.ions_native import main as ions_main


def _make_top_gro(tmp: Path):
    top = tmp / "a.top"
    gro = tmp / "a.gro"
    top.write_text("""[ defaults ]
1 1 yes 1 1
[ atomtypes ]
CA 1 0 A 0.1 0.2
[ moleculetype ]
Macro 1
[ atoms ]
1 CA 1 RES CA 1
[ bonds ]
[ angles ]
[ dihedrals ]
[ pairs ]
[ exclusions ]
[ system ]
X
[ molecules ]
Macro 1
""")
    gro.write_text("""h
1
    1RES     CA    1   0.000   0.000   0.000
1 1 1
""")
    return top, gro


def test_ions_with_explicit_params(tmp_path: Path):
    top, gro = _make_top_gro(tmp_path)
    out_top = tmp_path / "o.top"
    out_gro = tmp_path / "o.gro"
    rc = ions_main(["-f", str(top), "-g", str(gro), "-of", str(out_top), "-og", str(out_gro), "-ionnm", "K", "-ionn", "2", "-ionq", "1", "-ionm", "1", "-ionC12", "0.2", "-ionC6", "0.1"])
    assert rc == 0
    txt = out_top.read_text()
    assert "K 2" in txt
    assert "K\t1.0\t1.0\tA\t0.1\t0.2" in txt
    assert out_gro.read_text().splitlines()[1].strip() == "3"


def test_ions_with_template_file(tmp_path: Path):
    top, gro = _make_top_gro(tmp_path)
    tdir = tmp_path / "tmpl"
    tdir.mkdir()
    (tdir / "ions.def").write_text("NA 2.0 1.5 0.3 0.2\n")
    out_top = tmp_path / "t.top"
    out_gro = tmp_path / "t.gro"
    rc = ions_main(["-f", str(top), "-g", str(gro), "-of", str(out_top), "-og", str(out_gro), "-ionnm", "NA", "-ionn", "1", "-t", str(tdir)])
    assert rc == 0
    assert "NA 1" in out_top.read_text()
