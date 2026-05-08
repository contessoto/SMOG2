from pathlib import Path
from unittest.mock import patch

from smog3 import smogcheck_parity


def test_smogcheck_subset_runner_executes_with_mocked_backend(tmp_path: Path):
    report = tmp_path / "report.json"
    fake_report = {
        "ok": True,
        "cases": [{"case": 1, "status": "PASS", "feature_group": "topology parameter parity"}],
        "summary": {"PASS": 1},
        "feature_groups": {"topology parameter parity": {"PASS": 1}},
    }

    with patch("smog3.smogcheck_parity.run_campaign", return_value=fake_report):
        rc = smogcheck_parity.main([
            "--cases",
            "1",
            "--report-json",
            str(report),
        ])

    assert rc == 0
