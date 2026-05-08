from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .parity_direct import compare_existing_dirs
from .selected_parity import _case56_shadow_free_setup

ROOT = Path(__file__).resolve().parents[2]
TESTLIST = ROOT / "SMOG-CHECK" / "share" / "settings" / "smog2.testlist"

KNOWN_MODELS = {
    "AA",
    "CA",
    "AA-2cg",
    "AA-nb-cr2",
    "AA-match",
    "AA-CC1",
    "AA-CCD",
    "AA-BOND",
    "AA-DIHE",
    "AA-DIHE4",
    "AA-interactive",
    "CA-interactive",
}

MODEL_FLAGS = {
    ("AA", False): ["-AA"],
    ("AA", True): ["-AAgaussian"],
    ("CA", False): ["-CA"],
    ("CA", True): ["-CAgaussian"],
    ("AA-2cg", False): ["-AA2cg"],
    ("AA-CC1", False): ["-AACC1"],
    ("AA-CCD", False): ["-AACCD"],
    ("AA-BOND", False): ["-AABOND"],
    ("AA-DIHE", False): ["-AADIHE"],
    ("AA-DIHE4", False): ["-AADIHE4"],
    ("AA-match", False): ["-AAMATCH"],
}


@dataclass(frozen=True)
class SmogcheckCase:
    case_id: int
    stem: str
    modifiers: tuple[str, ...]
    model: str
    contact_model: str
    params: tuple[str, ...]
    raw: str

    @property
    def pdb(self) -> str:
        return f"SMOG-CHECK/share/PDB.files/{self.stem}.pdb"

    @property
    def opensmog(self) -> bool:
        return "OpenSMOG" in self.modifiers

    @property
    def freecoor(self) -> bool:
        return "freecoor" in self.modifiers

    @property
    def interactive(self) -> bool:
        return self.model.endswith("-interactive")

    @property
    def gaussian(self) -> bool:
        return "gaussian" in self.contact_model

    @property
    def user_contacts(self) -> bool:
        return "userC" in self.contact_model


def parse_testlist(path: Path = TESTLIST) -> list[SmogcheckCase]:
    cases: list[SmogcheckCase] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";") or ";" not in line:
            continue
        left, case_txt = line.rsplit(";", 1)
        try:
            case_id = int(case_txt.strip())
        except ValueError:
            continue
        tokens = left.split()
        if len(tokens) < 3:
            continue
        stem = tokens[0]
        model_index = next((idx for idx, token in enumerate(tokens[1:], start=1) if token in KNOWN_MODELS), None)
        if model_index is None or model_index + 1 >= len(tokens):
            continue
        cases.append(
            SmogcheckCase(
                case_id=case_id,
                stem=stem,
                modifiers=tuple(tokens[1:model_index]),
                model=tokens[model_index],
                contact_model=tokens[model_index + 1],
                params=tuple(tokens[model_index + 2 :]),
                raw=line,
            )
        )
    return cases


def feature_group(case: SmogcheckCase) -> str:
    if case.interactive or case.freecoor:
        return "freecoor/interactive"
    if case.opensmog:
        return "OpenSMOG XML"
    if case.model == "CA" or case.model == "CA-interactive":
        return "CA coarse-graining"
    if case.contact_model == "shadow-free":
        return "shadow-free/custom topology"
    if case.contact_model in {"shadow", "cutoff", "cutoff-gaussian", "shadow-match"}:
        return "template/map variants"
    if case.user_contacts:
        return "user contacts"
    if case.model not in {"AA", "CA"}:
        return "template/map variants"
    return "topology parameter parity"


def _rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _shell_join(args: list[str]) -> str:
    return " ".join(subprocess.list2cmdline([arg]) for arg in args)


def _model_flags(case: SmogcheckCase) -> list[str] | None:
    return MODEL_FLAGS.get((case.model, case.gaussian))


def _supported(case: SmogcheckCase) -> tuple[bool, str]:
    if case.interactive:
        return False, "interactive prompt workflow is not implemented natively"
    if case.freecoor:
        return False, "freecoor input normalization is not implemented natively"
    if case.model == "AA-nb-cr2":
        return False, "AA-nb-cr2 template is not implemented natively"
    if case.contact_model == "shadow-match":
        return False, "AA-match shadow-match template setup is not implemented natively"
    if case.contact_model in {"shadow", "cutoff", "cutoff-gaussian"}:
        return False, "parameterized non-default contact/template setup is not fully implemented natively"
    if case.contact_model == "shadow-free" and case.case_id != 56:
        return False, "shadow-free topology template parity is only partially implemented"
    if _model_flags(case) is None:
        return False, f"model {case.model} is not implemented natively"
    return True, ""


def _baseline_setup(case: SmogcheckCase) -> str:
    if case.case_id == 56:
        return _case56_shadow_free_setup()
    return ""


def _baseline_args(case: SmogcheckCase, outdir: Path) -> list[str] | None:
    flags = _model_flags(case)
    if case.case_id == 56:
        flags = ["-t", "temp.bifsif"]
    if flags is None:
        return None
    args = [
        "smog2",
        "-i",
        f"/workdir/{case.pdb}",
        "-SCMorig",
        "-keep4SCM",
        "-o",
        f"/workdir/{_rel(outdir / 'model.top')}",
        "-g",
        f"/workdir/{_rel(outdir / 'model.gro')}",
        "-n",
        f"/workdir/{_rel(outdir / 'model.ndx')}",
        "-s",
        f"/workdir/{_rel(outdir / 'model.contacts')}",
    ]
    if case.opensmog:
        args.extend(["-OpenSMOG", "-OpenSMOGxml", f"/workdir/{_rel(outdir / 'model.xml')}"])
    args.extend(flags)
    if case.user_contacts:
        args.extend(["-c", f"/workdir/SMOG-CHECK/share/PDB.files/{case.stem}.contacts"])
    return args


def _candidate_args(case: SmogcheckCase, outdir: Path) -> list[str] | None:
    flags = _model_flags(case)
    if flags is None:
        return None
    args = [
        "python3",
        "-c",
        "from smog3.smog2_native import main; import sys; raise SystemExit(main(sys.argv[1:]))",
        "-i",
        case.pdb,
        "-o",
        _rel(outdir / "model.top"),
        "-g",
        _rel(outdir / "model.gro"),
        "-n",
        _rel(outdir / "model.ndx"),
        "-s",
        _rel(outdir / "model.contacts"),
    ]
    if case.opensmog:
        args.extend(["-OpenSMOG", "-OpenSMOGxml", _rel(outdir / "model.xml")])
    args.extend(flags)
    if case.user_contacts:
        args.extend(["-c", f"SMOG-CHECK/share/PDB.files/{case.stem}.contacts"])
    if case.contact_model == "shadow-free":
        args.extend(["-contactMode", "shadow-free", "-contactParam", case.params[0] if case.params else "5.0"])
    return args


def _run_baseline(case: SmogcheckCase, outdir: Path, image: str) -> tuple[int, str, list[str]]:
    args = _baseline_args(case, outdir)
    if args is None:
        return 99, "unsupported baseline argument mapping", []
    setup = _baseline_setup(case)
    script = f"""
set -euo pipefail
cd /workdir
cd {_rel(outdir)}
{setup}
    {_shell_join(args)}
"""
    proc = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{ROOT}:/workdir", image, "bash", "-lc", script],
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr, args


def _run_candidate(case: SmogcheckCase, outdir: Path) -> tuple[int, str, list[str]]:
    args = _candidate_args(case, outdir)
    if args is None:
        return 99, "unsupported candidate argument mapping", []
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    env["SMOG3_USE_SCM_DEFAULTS"] = "1"
    proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout + proc.stderr, args


def _line_counts(directory: Path, include_xml: bool) -> dict[str, int | None]:
    names = ["model.top", "model.gro", "model.ndx", "model.contacts"] + (["model.xml"] if include_xml else [])
    return {name: (len((directory / name).read_text().splitlines()) if (directory / name).exists() else None) for name in names}


def _case_status(baseline_rc: int, candidate_rc: int, comparison_ok: bool) -> str:
    if baseline_rc != 0:
        return "BASELINE_ERROR"
    if candidate_rc != 0:
        return "CANDIDATE_ERROR"
    return "PASS" if comparison_ok else "FAIL"


def _select_cases(all_cases: list[SmogcheckCase], ns: argparse.Namespace) -> list[SmogcheckCase]:
    if ns.cases:
        wanted = {int(x) for x in ns.cases.split(",") if x.strip()}
        selected = [case for case in all_cases if case.case_id in wanted]
    elif ns.range:
        lo_txt, hi_txt = ns.range.split("-", 1)
        lo, hi = int(lo_txt), int(hi_txt)
        selected = [case for case in all_cases if lo <= case.case_id <= hi]
    elif ns.all:
        selected = list(all_cases)
    else:
        raise SystemExit("Choose --cases, --range, or --all")
    if ns.stop_after is not None:
        selected = selected[: ns.stop_after]
    return selected


def run_campaign(cases: list[SmogcheckCase], out_root: Path, image: str) -> dict:
    if shutil.which("docker") is None:
        return {"ok": False, "reason": "docker not available", "cases": []}
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True)
    reports_dir = out_root / "reports"
    reports_dir.mkdir()

    report: dict = {
        "ok": True,
        "image": image,
        "total_selected": len(cases),
        "cases": [],
        "summary": {},
        "feature_groups": {},
    }
    for case in cases:
        supported, reason = _supported(case)
        case_root = out_root / f"case{case.case_id}"
        baseline_dir = case_root / "baseline"
        candidate_dir = case_root / "candidate"
        baseline_dir.mkdir(parents=True)
        candidate_dir.mkdir(parents=True)
        group = feature_group(case)
        if not supported:
            entry = {
                "case": case.case_id,
                "pdb": case.stem,
                "model": case.model,
                "contact_model": case.contact_model,
                "feature_group": group,
                "status": "UNSUPPORTED",
                "reason": reason,
                "raw": case.raw,
            }
            report["cases"].append(entry)
            (reports_dir / f"case{case.case_id}.json").write_text(json.dumps(entry, indent=2))
            report["ok"] = False
            print(f"case {case.case_id}: UNSUPPORTED ({group})", flush=True)
            continue

        brc, bout, bargs = _run_baseline(case, baseline_dir, image)
        crc, cout, cargs = _run_candidate(case, candidate_dir)
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=case.opensmog)
        status = _case_status(brc, crc, comparison["ok"])
        entry = {
            "case": case.case_id,
            "pdb": case.stem,
            "model": case.model,
            "contact_model": case.contact_model,
            "feature_group": group,
            "status": status,
            "baseline_rc": brc,
            "candidate_rc": crc,
            "baseline_args": bargs,
            "candidate_args": cargs,
            "baseline_tail": bout[-2000:],
            "candidate_tail": cout[-2000:],
            "baseline_line_counts": _line_counts(baseline_dir, case.opensmog),
            "candidate_line_counts": _line_counts(candidate_dir, case.opensmog),
            "comparisons": comparison["comparisons"],
            "raw": case.raw,
        }
        report["cases"].append(entry)
        (reports_dir / f"case{case.case_id}.json").write_text(json.dumps(entry, indent=2))
        report["ok"] = bool(report["ok"] and status == "PASS")
        print(f"case {case.case_id}: {status} ({group})", flush=True)

    summary: dict[str, int] = {}
    groups: dict[str, dict[str, int]] = {}
    for entry in report["cases"]:
        status = entry["status"]
        summary[status] = summary.get(status, 0) + 1
        group = entry["feature_group"]
        groups.setdefault(group, {})
        groups[group][status] = groups[group].get(status, 0) + 1
    report["summary"] = summary
    report["feature_groups"] = groups
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--cases")
    group.add_argument("--range")
    group.add_argument("--all", action="store_true")
    parser.add_argument("--stop-after", type=int, default=None)
    parser.add_argument("--out-root", default="parity_runs/all")
    parser.add_argument("--report-json", default="parity_all_summary.json")
    parser.add_argument("--image", default="smogserver/smog2:stable")
    ns = parser.parse_args(argv)

    all_cases = parse_testlist()
    selected = _select_cases(all_cases, ns)
    report = run_campaign(selected, ROOT / ns.out_root, ns.image)
    report["discovered"] = {
        "total_cases": len(all_cases),
        "direct_smog2_cases": len(all_cases),
        "helper_tool_cases": 0,
    }
    (ROOT / ns.report_json).write_text(json.dumps(report, indent=2))

    print("case status feature")
    for entry in report.get("cases", []):
        print(f"{entry['case']:>4} {entry['status']:<15} {entry['feature_group']}")
    print("\nsummary")
    for status in ["PASS", "FAIL", "UNSUPPORTED", "BASELINE_ERROR", "CANDIDATE_ERROR"]:
        print(f"  {status:<15} {report.get('summary', {}).get(status, 0)}")
    print(f"\nreport written to {ns.report_json}")
    return 0 if report.get("ok") else 2


if __name__ == "__main__":
    raise SystemExit(main())
