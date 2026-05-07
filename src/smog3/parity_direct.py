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
    xml_args = ["-OpenSMOG", "-OpenSMOGxml", "model.xml"] if case_id == 94 else []
    userc = ["-userContacts", str(ROOT / "SMOG-CHECK" / "share" / "PDB.files" / "2ci2_v2.contacts")] if case_id == 50 else []
    args = [
        "python",
        "-c",
        "from smog3.smog2_native import main; import sys; raise SystemExit(main(sys.argv[1:]))",
        "-i",
        str(pdb),
        "-o",
        "model.top",
        "-g",
        "model.gro",
        "-n",
        "model.ndx",
        "-s",
        "model.contacts",
        *CASE_ARGS[case_id][:1],
        *userc,
        *xml_args,
    ]
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


def compare_existing_dirs(baseline_dir: Path, candidate_dir: Path, include_xml: bool = False) -> dict:
    files = ["model.top", "model.gro", "model.ndx", "model.contacts"] + (["model.xml"] if include_xml else [])
    comps = {f: _compare_file(baseline_dir / f, candidate_dir / f) for f in files}
    return {
        "baseline_dir": str(baseline_dir),
        "candidate_dir": str(candidate_dir),
        "comparisons": comps,
        "ok": all(v.get("match") for v in comps.values()),
    }


def run_cases(case_ids: list[int]) -> dict:
    if not _perl_ready():
        return {"skipped": True, "reason": "Perl dependencies missing (XML::Simple/XML::Validator::Schema/PDL)"}
    report = {"skipped": False, "cases": []}
    with tempfile.TemporaryDirectory(prefix="smog3-parity-") as td:
        root = Path(td)
        for cid in case_ids:
            bdir = root / f"baseline_{cid}"
            cdir = root / f"candidate_{cid}"
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
    p.add_argument("--compare-existing", action="store_true")
    p.add_argument("--baseline-dir", default=None)
    p.add_argument("--candidate-dir", default=None)
    p.add_argument("--include-xml", action="store_true")
    ns = p.parse_args(argv)

    if ns.compare_existing:
        if not ns.baseline_dir or not ns.candidate_dir:
            raise SystemExit("--compare-existing requires --baseline-dir and --candidate-dir")
        report = compare_existing_dirs(Path(ns.baseline_dir), Path(ns.candidate_dir), include_xml=ns.include_xml)
        Path(ns.report_json).write_text(json.dumps(report, indent=2))
        for f, comp in report["comparisons"].items():
            print(f"{f}: {'OK' if comp.get('match') else 'DIFF'}")
        return 0 if report.get("ok") else 2

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
