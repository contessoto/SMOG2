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
    assert "K\t1\t1\tA\t0.1\t0.2" in txt
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


def test_ions_with_template_named_ions_def_and_opensmog_xml(tmp_path: Path):
    top, gro = _make_top_gro(tmp_path)
    xml = tmp_path / "in.xml"
    xml.write_text(
        """<OpenSMOGforces>
 <nonbond>
  <nonbond_bytype>
   <expression expr="null"/>
   <parameter>null</parameter>
   <nonbond_param type1="NB_1" type2="NB_1" null="0"/>
  </nonbond_bytype>
 </nonbond>
</OpenSMOGforces>
"""
    )
    tdir = tmp_path / "tmpl"
    tdir.mkdir()
    (tdir / "custom.ions.def").write_text("MG 1 2 5.96046e-09 0\n")
    out_top = tmp_path / "t.top"
    out_gro = tmp_path / "t.gro"
    out_xml = tmp_path / "out.xml"

    rc = ions_main([
        "-f", str(top),
        "-g", str(gro),
        "-of", str(out_top),
        "-og", str(out_gro),
        "-OpenSMOG", str(xml),
        "-OpenSMOGout", str(out_xml),
        "-ionnm", "MG",
        "-ionn", "2",
        "-t", str(tdir),
    ])

    assert rc == 0
    top_text = out_top.read_text()
    assert top_text.index("[ moleculetype ]\nMG 1") < top_text.index("[ moleculetype ]\nMacro 1")
    assert out_gro.read_text().splitlines()[1].strip() == "3"
    xml_text = out_xml.read_text()
    assert 'type1="MG"' in xml_text or 'type2="MG"' in xml_text
    assert 'null="0.00000e+00"' not in xml_text
    assert 'null="0"' in xml_text
