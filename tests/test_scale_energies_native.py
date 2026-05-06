from pathlib import Path

from smog3.scale_energies_native import ScaleConfig, main as scale_main, scale_topology


def _write_inputs(tmp_path: Path):
    top_in = tmp_path / "in.top"
    ndx = tmp_path / "in.ndx"
    top_out = tmp_path / "out.top"

    top_in.write_text(
        """[ defaults ]
1 1 yes 1 1
[ pairs ]
1 2 1 0.5 0.1
2 3 1 0.5 0.1
[ dihedrals ]
1 2 3 4 1 180 2.0 1
1 2 3 5 1 180 2.0 1
[ exclusions ]
1 2
2 5
"""
    )
    ndx.write_text(
        """[ DG ]
1 2 3 4
[ C1 ]
1 2
[ C2 ]
2 3
"""
    )
    return top_in, ndx, top_out


def test_scale_topology_pairs_and_dihedrals(tmp_path: Path):
    top_in, ndx, top_out = _write_inputs(tmp_path)

    cfg = ScaleConfig(group_d="DG", group_c1="C1", group_c2="C2", rescale_c=2.0, rescale_d=3.0)
    scale_topology(top_in, ndx, top_out, cfg)

    txt = top_out.read_text()
    assert "1\t2\t1\t1.0\t0.2" in txt
    assert "2\t3\t1\t1.0\t0.2" in txt
    assert "1\t2\t3\t4\t1\t180\t6.0\t1" in txt
    assert "1 2 3 5 1 180 2.0 1" in txt


def test_scale_main_cli(tmp_path: Path):
    top_in, ndx, top_out = _write_inputs(tmp_path)
    rc = scale_main([
        "-f", str(top_in), "-n", str(ndx), "-of", str(top_out),
        "-rc", "0", "-rd", "2", "-grpD", "DG", "-grpC1", "C1", "-grpC2", "C2"
    ])
    assert rc == 0
    txt = top_out.read_text()
    assert "removed by smog_scale-energies" in txt
