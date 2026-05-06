import math
from pathlib import Path

from smog3.tablegen import main as tablegen_main


def _rows(path: Path):
    out = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            out.append([float(x) for x in line.split()])
    return out


def _assert_rows_close(a: Path, b: Path):
    ra = _rows(a)
    rb = _rows(b)
    assert len(ra) == len(rb)
    for xa, xb in zip(ra, rb):
        assert len(xa) == len(xb)
        for va, vb in zip(xa, xb):
            assert math.isclose(va, vb, rel_tol=1e-12, abs_tol=1e-12)


def test_tablegen_matches_repo_reference_kcal300(tmp_path: Path):
    out = tmp_path / "table.xvg"
    rc = tablegen_main(["-table", str(out), "-unit", "kCal", "-temp", "300", "-N", "6", "-M", "12", "-ic", "0", "-sd", "1.0", "-sc", "1.5", "-tl", "5.0"])
    assert rc == 0
    ref = Path("SMOG-CHECK/share/refs/table.kCal300.xvg")
    _assert_rows_close(ref, out)


def test_tablegen_matches_repo_reference_kj300(tmp_path: Path):
    out = tmp_path / "table.xvg"
    rc = tablegen_main(["-table", str(out), "-unit", "kJ", "-temp", "300", "-N", "6", "-M", "12", "-ic", "0", "-sd", "1.0", "-sc", "1.5", "-tl", "5.0"])
    assert rc == 0
    ref = Path("SMOG-CHECK/share/refs/table.kJ300.xvg")
    _assert_rows_close(ref, out)


def test_tablegen_matches_repo_reference_kj100(tmp_path: Path):
    out = tmp_path / "table.xvg"
    rc = tablegen_main(["-table", str(out), "-unit", "kJ", "-temp", "100", "-N", "6", "-M", "12", "-ic", "0", "-sd", "1.0", "-sc", "1.5", "-tl", "5.0"])
    assert rc == 0
    ref = Path("SMOG-CHECK/share/refs/table.kJ100.xvg")
    _assert_rows_close(ref, out)
