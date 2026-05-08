from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .parity_direct import compare_existing_dirs

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class SelectedCase:
    case_id: int
    stem: str
    baseline_args: tuple[str, ...]
    candidate_args: tuple[str, ...]
    include_xml: bool = False
    baseline_setup: str = ""

    @property
    def pdb(self) -> str:
        return f"SMOG-CHECK/share/PDB.files/{self.stem}.pdb"


def _case56_shadow_free_setup() -> str:
    # Mirrors the SMOG-CHECK case-56 template specialization for:
    # 1AKEapo_v2 AA shadow-free 5.0 1.4 1.2 1.0 2.0 1.0 1.0 1.0
    # 2.5 0.01 1.0 1.4 5.0 0.2 -1 1E-6 3E-9 1.0
    sigma_nm = 2.5 / 10.0
    epsilon = 0.01
    c12_nb1 = (sigma_nm**12) * epsilon
    template = "/opt/smog2/SMOG-CHECK/share/templates/SBM_AA"
    return f"""
mkdir -p temp.bifsif
sed "s/PARM_C_D/1.2/g;s/PARM_P_BB/1.0/g;s/PARM_P_SC/1.0/g;s/PARM_N_BB/1.0/g;s/PARM_N_SC/2.0/g;s/CUTDIST/5.0/g;s/SCM_R/1.4/g;s/SCM_BR/0.5/g;s/MINVERSION/2.4.5/g" {template}/AA-whitford09.shadow.free.sif > temp.bifsif/tmp.sif
sed "s/PARM_MASS/0.2/g;s/PARM_chargeNB/-1/g;s/PARM_C6_2/1E-6/g;s/PARM_C12_2/3E-9/g;s/PARM_C12/{c12_nb1:.12g}/g" {template}/AA-whitford09.free.nb > temp.bifsif/tmp.nb
cp {template}/AA-whitford09.free.bif temp.bifsif/tmp.bif
cp {template}/AA-whitford09.free.b temp.bifsif/tmp.b
cp {template}/extras temp.bifsif/test.extras
"""


SELECTED_CASES: dict[int, SelectedCase] = {
    1: SelectedCase(1, "1A01-AMP", ("-AA",), ("-AA",)),
    21: SelectedCase(21, "1A01-AMP", ("-AAgaussian",), ("-AAgaussian",)),
    41: SelectedCase(41, "1reschains_v2", ("-CA",), ("-CA",)),
    50: SelectedCase(
        50,
        "2ci2_v2",
        ("-AA", "-c", "/workdir/SMOG-CHECK/share/PDB.files/2ci2_v2.contacts"),
        ("-AA", "-c", "SMOG-CHECK/share/PDB.files/2ci2_v2.contacts"),
    ),
    56: SelectedCase(
        56,
        "1AKEapo_v2",
        ("-t", "temp.bifsif"),
        ("-AA", "-contactMode", "shadow-free", "-contactParam", "5.0"),
        baseline_setup=_case56_shadow_free_setup(),
    ),
    94: SelectedCase(
        94,
        "1AKEapo_v2",
        ("-AA", "-OpenSMOG", "-OpenSMOGxml", "/workdir/parity_runs/selected/case94/baseline/model.xml"),
        ("-AA", "-OpenSMOG", "-OpenSMOGxml", "parity_runs/selected/case94/candidate/model.xml"),
        include_xml=True,
    ),
}


def _rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _base_smog2_args(case: SelectedCase, outdir: Path) -> list[str]:
    return [
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
        *case.baseline_args,
    ]


def _candidate_args(case: SelectedCase, outdir: Path) -> list[str]:
    return [
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
        *case.candidate_args,
    ]


def _shell_join(args: list[str]) -> str:
    return " ".join(subprocess.list2cmdline([arg]) for arg in args)


def _run_baseline(case: SelectedCase, outdir: Path, image: str) -> tuple[int, str, list[str]]:
    args = _base_smog2_args(case, outdir)
    setup = case.baseline_setup
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


def _run_candidate(case: SelectedCase, outdir: Path) -> tuple[int, str, list[str]]:
    args = _candidate_args(case, outdir)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    env["SMOG3_USE_SCM_DEFAULTS"] = "1"
    proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout + proc.stderr, args


def _line_counts(directory: Path, include_xml: bool) -> dict[str, int | None]:
    files = ["model.top", "model.gro", "model.ndx", "model.contacts"] + (["model.xml"] if include_xml else [])
    out: dict[str, int | None] = {}
    for name in files:
        path = directory / name
        out[name] = len(path.read_text().splitlines()) if path.exists() else None
    return out


def _summarize_failures(compare: dict) -> dict[str, str]:
    out = {}
    for file_name, result in compare["comparisons"].items():
        if result.get("match"):
            continue
        reason = result.get("reason")
        if reason:
            out[file_name] = str(reason)
        else:
            diff = str(result.get("diff", ""))
            out[file_name] = diff.splitlines()[0] if diff else "diff"
    return out


def run_selected(case_ids: list[int], out_root: Path, image: str) -> dict:
    if shutil.which("docker") is None:
        return {"ok": False, "reason": "docker not available", "cases": []}

    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True)

    report = {"ok": True, "image": image, "cases": []}
    for case_id in case_ids:
        case = SELECTED_CASES[case_id]
        case_root = out_root / f"case{case_id}"
        baseline_dir = case_root / "baseline"
        candidate_dir = case_root / "candidate"
        baseline_dir.mkdir(parents=True)
        candidate_dir.mkdir(parents=True)

        brc, bout, bargs = _run_baseline(case, baseline_dir, image)
        crc, cout, cargs = _run_candidate(case, candidate_dir)
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=case.include_xml)
        case_ok = brc == 0 and crc == 0 and comparison["ok"]
        report["ok"] = bool(report["ok"] and case_ok)
        report["cases"].append(
            {
                "case": case_id,
                "pdb": case.stem,
                "baseline_rc": brc,
                "candidate_rc": crc,
                "baseline_args": bargs,
                "candidate_args": cargs,
                "baseline_tail": bout[-2000:],
                "candidate_tail": cout[-2000:],
                "baseline_line_counts": _line_counts(baseline_dir, case.include_xml),
                "candidate_line_counts": _line_counts(candidate_dir, case.include_xml),
                "comparisons": comparison["comparisons"],
                "failed_files": _summarize_failures(comparison),
                "ok": case_ok,
            }
        )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="1,21,41,50,56,94")
    parser.add_argument("--out-root", default="parity_runs/selected")
    parser.add_argument("--report-json", default="parity_selected.json")
    parser.add_argument("--image", default="smogserver/smog2:stable")
    ns = parser.parse_args(argv)

    case_ids = [int(x) for x in ns.cases.split(",") if x.strip()]
    unknown = [case_id for case_id in case_ids if case_id not in SELECTED_CASES]
    if unknown:
        raise SystemExit(f"Unsupported selected cases: {unknown}")

    report = run_selected(case_ids, ROOT / ns.out_root, ns.image)
    (ROOT / ns.report_json).write_text(json.dumps(report, indent=2))
    for case in report.get("cases", []):
        print(f"case {case['case']}: {'OK' if case['ok'] else 'DIFF'}")
        print(f"  baseline={case['baseline_rc']} candidate={case['candidate_rc']}")
        for file_name, result in case["comparisons"].items():
            print(f"  {file_name}: {'OK' if result.get('match') else 'DIFF'}")
    print(f"selected parity report written to {ns.report_json}")
    return 0 if report.get("ok") else 2


if __name__ == "__main__":
    raise SystemExit(main())
