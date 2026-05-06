from pathlib import Path

from smog3.extract_native import main as extract_main


def _make_inputs(tmp: Path):
    top = tmp / "in.top"
    gro = tmp / "in.gro"
    ndx = tmp / "in.ndx"
    top.write_text(
        """[ atoms ]
1 A 1 RES A1 1 0 1
2 A 1 RES A2 2 0 1
3 A 1 RES A3 3 0 1
4 A 1 RES A4 4 0 1
[ bonds ]
1 2 1
2 3 1
3 4 1
[ angles ]
1 2 3 1
2 3 4 1
[ dihedrals ]
1 2 3 4 1 180 1 1
[ pairs ]
1 3 1 0.5 0.1
2 4 1 0.5 0.1
[ exclusions ]
1 3
2 4
"""
    )
    gro.write_text(
        """header
4
    1RES     A1    1   0.000   0.000   0.000
    1RES     A2    2   0.100   0.000   0.000
    1RES     A3    3   1.000   0.000   0.000
    1RES     A4    4   2.000   0.000   0.000
1 1 1
"""
    )
    ndx.write_text("""[ SEL ]
1 2 3
""")
    return top, gro, ndx


def test_extract_with_ndx_group(tmp_path: Path):
    top, gro, ndx = _make_inputs(tmp_path)
    out_top = tmp_path / "out.top"
    out_gro = tmp_path / "out.gro"
    out_map = tmp_path / "out.map"

    rc = extract_main(["-f", str(top), "-g", str(gro), "-n", str(ndx), "-of", str(out_top), "-og", str(out_gro), "-om", str(out_map), "-group", "SEL"])
    assert rc == 0
    txt = out_top.read_text()
    assert "4" not in txt
    assert "1\t2\t1" in txt
    assert "2\t3\t1" in txt
    assert "3\t4\t1" not in txt
    assert out_gro.read_text().splitlines()[1].strip() == "3"


def test_extract_with_nondx_distance_selection(tmp_path: Path):
    top, gro, ndx = _make_inputs(tmp_path)
    out_top = tmp_path / "out2.top"
    out_gro = tmp_path / "out2.gro"
    out_map = tmp_path / "out2.map"

    rc = extract_main(["-f", str(top), "-g", str(gro), "-nondx", "-distfrom", "1", "-distval", "0.3", "-of", str(out_top), "-og", str(out_gro), "-om", str(out_map)])
    assert rc == 0
    assert out_gro.read_text().splitlines()[1].strip() == "2"

def test_extract_generates_restraints_and_xml_copy(tmp_path: Path):
    top, gro, ndx = _make_inputs(tmp_path)
    xml = tmp_path / "in.xml"
    xml.write_text("<OpenSMOG/>\n")
    out_top = tmp_path / "out3.top"
    out_gro = tmp_path / "out3.gro"
    out_map = tmp_path / "out3.map"
    out_xml = tmp_path / "out3.xml"
    rmap = tmp_path / "restrained.map"
    rc = extract_main([
        "-f", str(top), "-g", str(gro), "-n", str(ndx), "-group", "SEL",
        "-of", str(out_top), "-og", str(out_gro), "-om", str(out_map),
        "-OpenSMOG", str(xml), "-OpenSMOGout", str(out_xml), "-restraints", "100", "-orm", str(rmap)
    ])
    assert rc == 0
    assert out_xml.read_text() == "<OpenSMOG/>\n"
    assert rmap.exists()
