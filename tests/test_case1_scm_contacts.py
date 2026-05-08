from pathlib import Path
from types import SimpleNamespace

import pytest

from smog3 import smog2_native
from smog3.smog2_native import main as smog2_main


def _case1_pdb() -> Path:
    return Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files" / "1A01-AMP.pdb"


def test_case1_ndx_uses_scm_numeric_chain_groups(tmp_path: Path):
    base = tmp_path / "case1"
    rc = smog2_main([
        "-i", str(_case1_pdb()),
        "-AA",
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])

    assert rc == 0
    lines = base.with_suffix(".ndx").read_text().splitlines()
    assert [line for line in lines if line.startswith("[ ")] == ["[ 1 ]", "[ 2 ]", "[ 3 ]", "[ 4 ]"]
    assert lines[1:4] == ["1", "2", "3"]
    assert "[ System ]" not in lines
    assert "2702" in lines


def test_case1_scm_contact_generation_uses_java_not_perl(monkeypatch, tmp_path: Path):
    calls = []

    def fake_run(argv, capture_output, text):
        calls.append(argv)
        assert Path(argv[0]).name != "perl"
        out = Path(argv[argv.index("-o") + 1])
        out.write_text("1 2 3 4\n")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(smog2_native.shutil, "which", lambda name: f"/usr/bin/{name}" if name == "java" else None)
    monkeypatch.setattr(smog2_native.subprocess, "run", fake_run)

    base = tmp_path / "case1"
    rc = smog2_main([
        "-i", str(_case1_pdb()),
        "-AA",
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])

    assert rc == 0
    assert calls
    assert calls[0][0] == "/usr/bin/java"
    assert "-jar" in calls[0]
    assert calls[0][calls[0].index("-ch") + 1] == str(base.with_suffix(".ndx"))
    assert calls[0][calls[0].index("-t") + 1] == str(base.with_suffix(".top").with_name("case1.top4SCM.top"))
    for flag, value in [("-m", "shadow"), ("-c", "6"), ("-s", "1"), ("-br", "0.5"), ("-pd", "3")]:
        assert calls[0][calls[0].index(flag) + 1] == value
    assert "-bif" in calls[0]
    assert "--smog2output" in calls[0]
    assert "--showProgress" in calls[0]
    assert "--default" not in calls[0]
    assert base.with_suffix(".contacts").read_text() == "1 2 3 4\n"
    top4scm = base.with_suffix(".top").with_name("case1.top4SCM.top")
    assert top4scm.exists()
    top_text = top4scm.read_text()
    assert "[ bonds ]" in top_text
    assert "[ angles ]" in top_text
    assert "[ dihedrals ]" in top_text


def test_case1_contact_generation_nonempty_when_scm_jar_available(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    if not (root / "src" / "tools" / "SCM.jar").exists():
        pytest.skip("SCM.jar is not available")
    if smog2_native.shutil.which("java") is None:
        pytest.skip("java is not available")

    base = tmp_path / "case1"
    rc = smog2_main([
        "-i", str(_case1_pdb()),
        "-AA",
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])

    assert rc == 0
    contact_lines = [line for line in base.with_suffix(".contacts").read_text().splitlines() if line.strip()]
    assert contact_lines
    assert len(contact_lines) > 3306
    assert "1 53 4 2698" in contact_lines
    gro_lines = base.with_suffix(".gro").read_text().splitlines()
    assert gro_lines[26] == "    4LYS     CA   25  -0.914  -1.113   6.559"

    final_counts = {}
    current = None
    for raw in base.with_suffix(".top").read_text().splitlines():
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            final_counts[current] = 0
        elif current and s and not s.startswith(";"):
            final_counts[current] += 1
    assert final_counts["atoms"] == 2702
    assert final_counts["bonds"] == 2766
    assert final_counts["angles"] == 3737
    assert final_counts["dihedrals"] == 8260
    assert final_counts["pairs"] == len(contact_lines)
    assert final_counts["exclusions"] == len(contact_lines)
    pair_lines = []
    current = None
    for raw in base.with_suffix(".top").read_text().splitlines():
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            continue
        if current == "pairs" and s and not s.startswith(";"):
            pair_lines.append(s)
    assert pair_lines[0] == "1\t1195\t1\t 3.944505603e-02 7.999533051e-04"
    dihedral_lines = []
    current = None
    for raw in base.with_suffix(".top").read_text().splitlines():
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            continue
        if current == "dihedrals" and s and not s.startswith(";"):
            dihedral_lines.append(s)
    assert dihedral_lines[0] == "1\t2\t3\t4\t1\t 5.559814536e+02 2.303921569e-01 1"
    assert dihedral_lines[1] == "1\t2\t3\t4\t1\t 1.307944361e+03 1.151960784e-01 3"
    assert "200\t201\t203\t205\t2\t -7.620696300e-02 4.000000000e+01" in dihedral_lines

    top4scm = base.with_name("case1.top4SCM.top")
    counts = {}
    current = None
    for raw in top4scm.read_text().splitlines():
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            counts[current] = 0
        elif current and s and not s.startswith(";"):
            counts[current] += 1
    assert counts["atoms"] == 2702
    assert counts["bonds"] == 2766
    assert counts["angles"] == 3737
    assert counts["dihedrals"] == 8260
    gro4scm = base.with_name("case1.gro4SCM.gro")
    assert gro4scm.exists()
    lines = gro4scm.read_text().splitlines()
    assert lines[0] == "Temp Gro file with PDB precision for SCM calculations."
    assert int(lines[1]) == 2702
    assert lines[2].startswith("    1VAL      N    1")
    assert ".7853" in lines[2]
    assert lines[-1].split() == ["8.3908", "7.01", "9.5239"]


def test_case1_parity_script_cleans_outputs_and_prints_diagnostics():
    script = (Path(__file__).resolve().parents[1] / "scripts" / "run_case1_two_stage_parity.sh").read_text()
    assert 'rm -rf "$BASE_DIR" "$CAND_DIR"' in script
    assert 'rm -f "$REPORT_JSON"' in script
    assert 'mkdir -p "$BASE_DIR" "$CAND_DIR"' in script
    assert "ERROR: baseline output missing" in script
    assert "-keep4SCM" in script
    assert 'wc -l "$BASE_DIR/model.contacts" "$CAND_DIR/model.contacts"' in script
    assert 'exit "$COMPARE_RC"' in script


def test_selected_parity_script_uses_official_docker_and_selected_cases():
    script = (Path(__file__).resolve().parents[1] / "scripts" / "run_selected_two_stage_parity.sh").read_text()
    assert "smog3.selected_parity" in script
    assert '--cases "1,21,41,50,56,94"' in script
    assert "parity_selected.json" in script
