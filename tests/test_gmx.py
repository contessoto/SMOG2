from pathlib import Path

from smog3.gmx import parse_ndx, parse_top_sections, write_top_sections


def test_parse_ndx_groups(tmp_path: Path):
    ndx = tmp_path / "a.ndx"
    ndx.write_text("""[ Group1 ]\n1 2 3\n4\n\n[ Group2 ]\n7 8\n""")
    g = parse_ndx(ndx)
    assert g["Group1"] == [1, 2, 3, 4]
    assert g["Group2"] == [7, 8]


def test_parse_and_write_top_sections_roundtrip(tmp_path: Path):
    top = tmp_path / "a.top"
    text = "; preamble\n#include \"x.itp\"\n[ atoms ]\n1 A\n[ bonds ]\n1 2\n"
    top.write_text(text)
    secs = parse_top_sections(top)
    assert secs[1].name == "atoms"
    assert secs[2].name == "bonds"

    out = tmp_path / "b.top"
    write_top_sections(out, secs)
    assert out.read_text() == text
