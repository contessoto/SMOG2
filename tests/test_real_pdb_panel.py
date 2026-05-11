from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HELPER_PATH = ROOT / "scripts" / "compare_real_pdb_panel.py"
spec = importlib.util.spec_from_file_location("compare_real_pdb_panel", HELPER_PATH)
real_panel = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(real_panel)


def test_topology_section_counts_and_contact_counts(tmp_path: Path):
    top = tmp_path / "model.top"
    top.write_text(
        "[ atoms ]\n"
        "1 A 1 A A 1 0 1\n"
        "2 B 1 A B 1 0 1\n"
        "\n"
        "[ bonds ]\n"
        "; comment\n"
        "1 2 1 0.1 100\n",
        encoding="utf-8",
    )
    contacts = tmp_path / "model.contacts"
    contacts.write_text("; comment\n1 2 1 5\n\n2 3 1 6\n", encoding="utf-8")

    assert real_panel._topology_section_counts(top) == {"atoms": 2, "bonds": 1}
    assert real_panel._count_contacts(contacts) == 2


def test_gro_coordinate_delta_reports_max_difference(tmp_path: Path):
    baseline = tmp_path / "baseline.gro"
    candidate = tmp_path / "candidate.gro"
    baseline.write_text(
        "test\n"
        "2\n"
        "    1ALA     N    1   0.100   0.200   0.300\n"
        "    1ALA    CA    2   0.400   0.500   0.600\n"
        "   1.0 1.0 1.0\n",
        encoding="utf-8",
    )
    candidate.write_text(
        "test\n"
        "2\n"
        "    1ALA     N    1   0.100   0.201   0.300\n"
        "    1ALA    CA    2   0.400   0.500   0.602\n"
        "   1.0 1.0 1.0\n",
        encoding="utf-8",
    )

    delta = real_panel._gro_coordinate_delta(baseline, candidate)
    assert delta["available"] is True
    assert delta["baseline_atoms"] == 2
    assert delta["candidate_atoms"] == 2
    assert delta["differing_coordinate_fields"] == 2
    assert abs(delta["max_abs_coordinate_delta"] - 0.002) < 1e-12


def test_summary_counts_and_markdown():
    reports = [
        {"case": "a", "status": "PASS"},
        {"case": "b", "status": "DIFF"},
        {"case": "c", "status": "SMOG3_ERROR"},
    ]
    counts = real_panel._status_counts(reports)
    assert counts["PASS"] == 1
    assert counts["DIFF"] == 1
    assert counts["SMOG3_ERROR"] == 1

    markdown = real_panel._render_summary_markdown(
        {
            "smog3_version": "0.1.0a1",
            "docker_image": "smogserver/smog2:stable",
            "smog2_docker_version": "SMOG 2",
            "smog3_perl_invocations": 0,
            "counts": counts,
            "cases": [{"case": "a", "status": "PASS", "input_source": "fallback", "report": "reports/a.json"}],
        }
    )
    assert "Real PDB Panel Validation" in markdown
    assert "| a | PASS | fallback | `reports/a.json` |" in markdown


def test_download_error_status(tmp_path: Path):
    report = real_panel.build_case_report(
        case_name="bad_download",
        input_pdb=tmp_path / "missing.pdb",
        input_source="download-error",
        baseline_dir=tmp_path / "baseline",
        candidate_dir=tmp_path / "candidate",
        baseline_rc=99,
        candidate_rc=99,
        baseline_log=tmp_path / "baseline.log",
        candidate_log=tmp_path / "candidate.log",
        include_xml=False,
    )
    assert report["status"] == "DOWNLOAD_ERROR"
