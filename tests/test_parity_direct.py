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
