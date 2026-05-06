from pathlib import Path

from smog3.modifyxml_native import main as mod_main


def test_modifyxml_scale_param(tmp_path: Path):
    inp = tmp_path / "in.xml"
    out = tmp_path / "out.xml"
    inp.write_text('<interaction group1="A" group2="B" epsilon="2.0"/>\n')
    rc = mod_main(["-OpenSMOG", str(inp), "-OpenSMOGout", str(out), "-modgrp1", "A", "-modgrp2", "B", "-modparam", "epsilon", "-modby", "2"])
    assert rc == 0
    assert 'epsilon="4.0"' in out.read_text()


def test_modifyxml_remove_group(tmp_path: Path):
    inp = tmp_path / "in.xml"
    out = tmp_path / "out.xml"
    inp.write_text('<interaction group1="A" epsilon="2.0"/>\n<interaction group1="C" epsilon="1.0"/>\n')
    rc = mod_main(["-OpenSMOG", str(inp), "-OpenSMOGout", str(out), "-modgrp1", "A", "-remove"])
    assert rc == 0
    txt = out.read_text()
    assert 'group1="A"' not in txt
    assert 'group1="C"' in txt
