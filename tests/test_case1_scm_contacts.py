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
    assert base.with_suffix(".contacts").read_text() == "1 2 3 4\n"


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
    assert base.with_name("case1.gro4SCM.gro").exists()


def test_case1_parity_script_cleans_outputs_and_prints_diagnostics():
    script = (Path(__file__).resolve().parents[1] / "scripts" / "run_case1_two_stage_parity.sh").read_text()
    assert 'rm -rf "$BASE_DIR" "$CAND_DIR"' in script
    assert 'rm -f "$REPORT_JSON"' in script
    assert 'mkdir -p "$BASE_DIR" "$CAND_DIR"' in script
    assert "ERROR: baseline output missing" in script
    assert 'wc -l "$BASE_DIR/model.contacts" "$CAND_DIR/model.contacts"' in script
    assert 'exit "$COMPARE_RC"' in script
