from pathlib import Path
from unittest.mock import patch

from smog3 import smogcheck_parity


def test_smogcheck_subset_runner_executes_with_mocked_backend(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    report = tmp_path / "report.json"

    with patch("smog3.smogcheck_parity._run_smogcheck", return_value=0):
        rc = smogcheck_parity.main([
            "--repo-root",
            str(repo_root),
            "--start",
            "1",
            "--end",
            "1",
            "--report-json",
            str(report),
        ])

    assert rc == 0
