import math
from pathlib import Path

from smog3.tablegen import main as tablegen_main


def _rows(path: Path):
    rows = []
    for ln in path.read_text().splitlines():
        if ln and not ln.startswith("#"):
            rows.append([float(x) for x in ln.split()])
    return rows


def test_tablegen_default_matches_distributed_reference(tmp_path: Path):
    out = tmp_path / "table.xvg"
    rc = tablegen_main(["-table", str(out)])
    assert rc == 0

    ref = Path("SMOG-CHECK/share/refs/table.kCal300.xvg")
    ro = _rows(out)
    rr = _rows(ref)
    assert len(ro) == len(rr)
    for a, b in zip(ro, rr):
        assert len(a) == len(b)
        for va, vb in zip(a, b):
            assert math.isclose(va, vb, rel_tol=1e-12, abs_tol=1e-12)
