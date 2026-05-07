from __future__ import annotations

import argparse
import difflib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CASE_ARGS = {
    1: ["-AA"],
    21: ["-AAgaussian"],
    41: ["-CA"],
    50: ["-AA", "-c", str(ROOT / "SMOG-CHECK" / "share" / "PDB.files" / "2ci2_v2.contacts")],
    56: ["-AA"],
    94: ["-AA", "-OpenSMOG", "-OpenSMOGxml", "model.xml"],
}
CASE_PDB = {
    1: "1A01-AMP.pdb",
    21: "1A01-AMP.pdb",
    41: "1reschains_v2.pdb",
    50: "2ci2_v2.pdb",
    56: "1AKEapo_v2.pdb",
    94: "1AKEapo_v2.pdb",
}


def _perl_ready() -> bool:
    probe = subprocess.run(["perl", "-MXML::Simple", "-MXML::Validator::Schema", "-MPDL", "-e", "print 'ok'"], capture_output=True)
    return probe.returncode == 0


def _run_baseline(case_id: int, outdir: Path) -> tuple[int, str]:
    pdb = ROOT / "SMOG-CHECK" / "share" / "PDB.files" / CASE_PDB[case_id]
    args = ["perl", str(ROOT / "src" / "smogv2"), "-i", str(pdb), "-o", "model.top", "-g", "model.gro", "-n", "model.ndx", "-s", "model.contacts", *CASE_ARGS[case_id]]
    p = subprocess.run(args, cwd=outdir, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def _run_candidate(case_id: int, outdir: Path) -> tuple[int, str]:
    pdb = ROOT / "SMOG-CHECK" / "share" / "PDB.files" / CASE_PDB[case_id]
    xml_args = []
    if case_id == 94:
        xml_args = ["-OpenSMOG", "-OpenSMOGxml", "model.xml"]
    userc = []
    if case_id == 50:
        userc = ["-userContacts", str(ROOT / "SMOG-CHECK" / "share" / "PDB.files" / "2ci2_v2.contacts")]
    args = ["python", "-c", "from smog3.smog2_native import main; import sys; raise SystemExit(main(sys.argv[1:]))", "-i", str(pdb), "-o", "model.top", "-g", "model.gro", "-n", "model.ndx", "-s", "model.contacts", *CASE_ARGS[case_id][:1], *userc, *xml_args]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    p = subprocess.run(args, cwd=outdir, capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def _compare_file(a: Path, b: Path) -> dict:
    if not a.exists() or not b.exists():
        return {"match": False, "reason": "missing file"}
    ta, tb = a.read_text().splitlines(), b.read_text().splitlines()
    if ta == tb:
        return {"match": True}
    diff = "\n".join(difflib.unified_diff(ta, tb, fromfile=str(a), tofile=str(b), n=2))
    return {"match": False, "diff": diff[:4000]}


def run_cases(case_ids: list[int]) -> dict:
    if not _perl_ready():
        return {"skipped": True, "reason": "Perl dependencies missing (XML::Simple/XML::Validator::Schema/PDL)"}
    report = {"skipped": False, "cases": []}
    with tempfile.TemporaryDirectory(prefix="smog3-parity-") as td:
        root = Path(td)
        for cid in case_ids:
            bdir = root / f"baseline_{cid}"; cdir = root / f"candidate_{cid}"
            bdir.mkdir(); cdir.mkdir()
            brc, bout = _run_baseline(cid, bdir)
            crc, cout = _run_candidate(cid, cdir)
            files = ["model.top", "model.gro", "model.ndx", "model.contacts"] + (["model.xml"] if cid == 94 else [])
            comps = {f: _compare_file(bdir / f, cdir / f) for f in files}
            report["cases"].append({"case": cid, "baseline_rc": brc, "candidate_rc": crc, "baseline_out": bout[-500:], "candidate_out": cout[-500:], "comparisons": comps})
    report["ok"] = all(c["baseline_rc"] == 0 and c["candidate_rc"] == 0 and all(v.get("match") for v in c["comparisons"].values()) for c in report["cases"])
    return report


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cases", default="1,21,41,50,56,94")
    p.add_argument("--report-json", default="parity_report.json")
    ns = p.parse_args(argv)
    cases = [int(x) for x in ns.cases.split(",") if x.strip()]
    report = run_cases(cases)
    Path(ns.report_json).write_text(json.dumps(report, indent=2))
    if report.get("skipped"):
        print(report["reason"])
        return 4
    for c in report["cases"]:
        print(f"case {c['case']}: baseline={c['baseline_rc']} candidate={c['candidate_rc']}")
        for f, comp in c["comparisons"].items():
            print(f"  {f}: {'OK' if comp.get('match') else 'DIFF'}")
    return 0 if report.get("ok") else 2


if __name__ == "__main__":
    raise SystemExit(main())
