from pathlib import Path

from smog3 import parity_direct


def test_parity_direct_skips_cleanly_without_perl(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(parity_direct, "_perl_ready", lambda: False)
    rc = parity_direct.main(["--cases", "1", "--report-json", str(tmp_path / "r.json")])
    assert rc == 4
    txt = (tmp_path / "r.json").read_text()
    assert "Perl dependencies missing" in txt


def test_parity_report_structure_when_run_cases_skipped(monkeypatch):
    monkeypatch.setattr(parity_direct, "_perl_ready", lambda: False)
    report = parity_direct.run_cases([1, 21])
    assert report["skipped"] is True
    assert "Perl dependencies missing" in report["reason"]


def test_compare_existing_dirs(tmp_path: Path):
    b = tmp_path / "b"; c = tmp_path / "c"
    b.mkdir(); c.mkdir()
    for fn in ["model.top", "model.gro", "model.ndx", "model.contacts"]:
        (b / fn).write_text("same\n")
        (c / fn).write_text("same\n")
    report = parity_direct.compare_existing_dirs(b, c)
    assert report["ok"] is True
    rc = parity_direct.main(["--compare-existing", "--baseline-dir", str(b), "--candidate-dir", str(c), "--report-json", str(tmp_path / "rep.json")])
    assert rc == 0


def test_compare_existing_reports_missing_candidate_files(tmp_path: Path):
    b = tmp_path / "b"; c = tmp_path / "c"
    b.mkdir(); c.mkdir()
    for fn in ["model.top", "model.gro", "model.ndx", "model.contacts"]:
        (b / fn).write_text("baseline\n")
    report = parity_direct.compare_existing_dirs(b, c)
    assert report["ok"] is False
    for fn in ["model.top", "model.gro", "model.ndx", "model.contacts"]:
        assert report["comparisons"][fn]["reason"] == "missing file"
