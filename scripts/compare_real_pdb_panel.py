#!/usr/bin/env python3
from __future__ import annotations

import argparse
import difflib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from smog3.parity_direct import compare_existing_dirs


EXPECTED_OUTPUTS = ("model.top", "model.gro", "model.ndx", "model.contacts")
COMPARATOR_POLICY = (
    "Exact file comparison, with the existing documented parity policy for harmless "
    "SMOG-generated topology/XML metadata and tiny topology floating-point print noise."
)


def _read_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def _first_diff(a: Path, b: Path, *, limit: int = 80) -> str:
    if not a.exists() or not b.exists():
        return ""
    diff = difflib.unified_diff(
        _read_lines(a),
        _read_lines(b),
        fromfile=str(a),
        tofile=str(b),
        n=3,
        lineterm="",
    )
    return "\n".join(list(diff)[:limit])


def _gro_atom_count(path: Path) -> int | None:
    lines = _read_lines(path)
    if len(lines) < 2:
        return None
    try:
        return int(lines[1].strip())
    except ValueError:
        return None


def _gro_coords(line: str) -> tuple[float, float, float] | None:
    if len(line) >= 44:
        try:
            return (float(line[20:28]), float(line[28:36]), float(line[36:44]))
        except ValueError:
            pass
    parts = line.split()
    if len(parts) >= 3:
        try:
            return (float(parts[-3]), float(parts[-2]), float(parts[-1]))
        except ValueError:
            return None
    return None


def _gro_coordinate_delta(a: Path, b: Path) -> dict[str, Any]:
    a_count = _gro_atom_count(a)
    b_count = _gro_atom_count(b)
    if a_count is None or b_count is None or not a.exists() or not b.exists():
        return {"available": False, "reason": "missing or malformed GRO atom count"}
    a_lines = _read_lines(a)[2 : 2 + a_count]
    b_lines = _read_lines(b)[2 : 2 + b_count]
    compared = min(len(a_lines), len(b_lines))
    max_delta = 0.0
    differing_fields = 0
    for left, right in zip(a_lines, b_lines):
        left_coords = _gro_coords(left)
        right_coords = _gro_coords(right)
        if left_coords is None or right_coords is None:
            continue
        for x, y in zip(left_coords, right_coords):
            delta = abs(x - y)
            max_delta = max(max_delta, delta)
            if delta:
                differing_fields += 1
    return {
        "available": True,
        "baseline_atoms": a_count,
        "candidate_atoms": b_count,
        "atom_lines_compared": compared,
        "coordinate_fields_compared": compared * 3,
        "differing_coordinate_fields": differing_fields,
        "max_abs_coordinate_delta": max_delta,
    }


def _count_contacts(path: Path) -> int | None:
    if not path.exists():
        return None
    return sum(1 for line in _read_lines(path) if line.strip() and not line.lstrip().startswith(";"))


def _topology_section_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    current: str | None = None
    for line in _read_lines(path):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current = stripped.strip("[]").strip()
            counts.setdefault(current, 0)
            continue
        if current and stripped and not stripped.startswith(";"):
            counts[current] = counts.get(current, 0) + 1
    return counts


def _file_diagnostics(baseline_dir: Path, candidate_dir: Path) -> dict[str, Any]:
    diagnostics: dict[str, Any] = {}
    for name in EXPECTED_OUTPUTS + ("model.xml",):
        left = baseline_dir / name
        right = candidate_dir / name
        if not left.exists() and not right.exists():
            continue
        item: dict[str, Any] = {
            "baseline_exists": left.exists(),
            "candidate_exists": right.exists(),
            "baseline_lines": len(_read_lines(left)) if left.exists() else None,
            "candidate_lines": len(_read_lines(right)) if right.exists() else None,
        }
        if left.exists() and right.exists() and left.read_bytes() != right.read_bytes():
            item["first_diff"] = _first_diff(left, right)
        if name == "model.gro":
            item["gro_coordinate_delta"] = _gro_coordinate_delta(left, right)
        if name == "model.contacts":
            item["baseline_contact_lines"] = _count_contacts(left)
            item["candidate_contact_lines"] = _count_contacts(right)
        if name == "model.top":
            item["baseline_section_counts"] = _topology_section_counts(left)
            item["candidate_section_counts"] = _topology_section_counts(right)
        diagnostics[name] = item
    return diagnostics


def build_case_report(
    *,
    case_name: str,
    input_pdb: Path,
    input_source: str,
    baseline_dir: Path,
    candidate_dir: Path,
    baseline_rc: int,
    candidate_rc: int,
    baseline_log: Path,
    candidate_log: Path,
    include_xml: bool,
) -> dict[str, Any]:
    if input_source == "download-error":
        status = "DOWNLOAD_ERROR"
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=include_xml)
    elif baseline_rc != 0:
        status = "SMOG2_ERROR"
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=include_xml)
    elif candidate_rc != 0:
        status = "SMOG3_ERROR"
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=include_xml)
    else:
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=include_xml)
        status = "PASS" if comparison.get("ok") else "DIFF"
    return {
        "case": case_name,
        "status": status,
        "input_pdb": str(input_pdb),
        "input_source": input_source,
        "baseline_dir": str(baseline_dir),
        "candidate_dir": str(candidate_dir),
        "baseline_returncode": baseline_rc,
        "candidate_returncode": candidate_rc,
        "baseline_log": str(baseline_log),
        "candidate_log": str(candidate_log),
        "comparator_policy": COMPARATOR_POLICY,
        "comparison": comparison,
        "diagnostics": _file_diagnostics(baseline_dir, candidate_dir),
    }


def _load_reports(reports_dir: Path) -> list[dict[str, Any]]:
    return [json.loads(path.read_text(encoding="utf-8")) for path in sorted(reports_dir.glob("*.json"))]


def _status_counts(reports: list[dict[str, Any]]) -> dict[str, int]:
    counter = Counter(report.get("status", "UNKNOWN") for report in reports)
    return {key: counter.get(key, 0) for key in ("PASS", "DIFF", "SMOG2_ERROR", "SMOG3_ERROR", "DOWNLOAD_ERROR", "UNSUPPORTED")}


def _render_summary_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Real PDB Panel Validation",
        "",
        f"- SMOG3 version: `{summary.get('smog3_version', 'unknown')}`",
        f"- SMOG2 Docker image: `{summary.get('docker_image', 'unknown')}`",
        f"- SMOG2 Docker version: `{summary.get('smog2_docker_version', 'unknown')}`",
        f"- SMOG3 Perl invocations: `{summary.get('smog3_perl_invocations', 'unknown')}`",
        f"- Comparator policy: {COMPARATOR_POLICY}",
        "",
        "| Status | Count |",
        "| --- | ---: |",
    ]
    for status, count in summary["counts"].items():
        lines.append(f"| {status} | {count} |")
    lines.extend(
        [
            "",
            "| Case | Status | Input Source | Report |",
            "| --- | --- | --- | --- |",
        ]
    )
    for report in summary["cases"]:
        lines.append(
            f"| {report['case']} | {report['status']} | {report.get('input_source', '')} | "
            f"`{report['report']}` |"
        )
    lines.append("")
    return "\n".join(lines)


def compare_main(args: argparse.Namespace) -> int:
    baseline_dir = Path(args.baseline_dir)
    candidate_dir = Path(args.candidate_dir)
    include_xml = args.include_xml or (baseline_dir / "model.xml").exists() or (candidate_dir / "model.xml").exists()
    report = build_case_report(
        case_name=args.case,
        input_pdb=Path(args.input_pdb),
        input_source=args.input_source,
        baseline_dir=baseline_dir,
        candidate_dir=candidate_dir,
        baseline_rc=args.baseline_rc,
        candidate_rc=args.candidate_rc,
        baseline_log=Path(args.baseline_log),
        candidate_log=Path(args.candidate_log),
        include_xml=include_xml,
    )
    Path(args.report_json).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"{args.case}: {report['status']}")
    return 0 if report["status"] == "PASS" else 1


def summarize_main(args: argparse.Namespace) -> int:
    reports_dir = Path(args.reports_dir)
    reports = _load_reports(reports_dir)
    cases = [
        {
            "case": report["case"],
            "status": report["status"],
            "input_source": report.get("input_source"),
            "report": str(reports_dir / f"{report['case']}.json"),
        }
        for report in reports
    ]
    summary = {
        "smog3_version": args.smog3_version,
        "smog2_docker_version": args.smog2_version,
        "docker_image": args.docker_image,
        "case_count": len(reports),
        "counts": _status_counts(reports),
        "smog3_perl_invocations": args.smog3_perl_invocations,
        "smog3_perl_log": args.smog3_perl_log,
        "comparator_policy": COMPARATOR_POLICY,
        "cases": cases,
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    Path(args.summary_md).write_text(_render_summary_markdown(summary), encoding="utf-8")
    print(f"summary written to {args.summary_json} and {args.summary_md}")
    return 0 if summary["counts"].get("PASS", 0) == len(reports) else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare and summarize real-PDB SMOG2/SMOG3 validation outputs.")
    sub = parser.add_subparsers(dest="command", required=True)

    compare = sub.add_parser("compare")
    compare.add_argument("--case", required=True)
    compare.add_argument("--input-pdb", required=True)
    compare.add_argument("--input-source", required=True)
    compare.add_argument("--baseline-dir", required=True)
    compare.add_argument("--candidate-dir", required=True)
    compare.add_argument("--baseline-rc", type=int, required=True)
    compare.add_argument("--candidate-rc", type=int, required=True)
    compare.add_argument("--baseline-log", required=True)
    compare.add_argument("--candidate-log", required=True)
    compare.add_argument("--report-json", required=True)
    compare.add_argument("--include-xml", action="store_true")
    compare.set_defaults(func=compare_main)

    summarize = sub.add_parser("summarize")
    summarize.add_argument("--reports-dir", required=True)
    summarize.add_argument("--summary-json", required=True)
    summarize.add_argument("--summary-md", required=True)
    summarize.add_argument("--smog3-version", required=True)
    summarize.add_argument("--smog2-version", required=True)
    summarize.add_argument("--docker-image", required=True)
    summarize.add_argument("--smog3-perl-invocations", type=int, required=True)
    summarize.add_argument("--smog3-perl-log", required=True)
    summarize.set_defaults(func=summarize_main)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
