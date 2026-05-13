from __future__ import annotations

import argparse
import difflib
import itertools
import json
import math
import os
import re
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
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    p = subprocess.run(args, cwd=outdir, capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def _drop_top_header_metadata(lines: list[str]) -> list[str]:
    for idx, line in enumerate(lines):
        if line.strip().startswith("["):
            return lines[idx:]
    return lines


def _drop_xml_header_metadata(lines: list[str]) -> list[str]:
    for idx, line in enumerate(lines):
        if line.strip().startswith("<OpenSMOGforces"):
            return lines[idx:]
    return lines


_TOP_FLOAT_SECTIONS = {"atomtypes", "bonds", "angles", "dihedrals", "pairs"}
_TOP_FLOAT_ABS_TOL = 5e-10
_TOP_FLOAT_REL_TOL = 5e-9


def _top_section(line: str) -> str | None:
    m = re.match(r"\s*\[\s*([^\]]+?)\s*\]", line)
    return m.group(1).strip().lower() if m else None


def _float_like(token: str) -> bool:
    return any(ch in token for ch in ".eE") and _as_float(token) is not None


def _as_float(token: str) -> float | None:
    try:
        value = float(token)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _split_top_data_comment(line: str) -> tuple[str, str]:
    idx = line.find(";")
    if idx < 0:
        return line, ""
    return line[:idx], line[idx:]


def _top_tokens_with_layout(data: str) -> tuple[list[tuple[str, str]], str]:
    out: list[tuple[str, str]] = []
    pos = 0
    for match in re.finditer(r"\S+", data):
        out.append((data[pos : match.start()], match.group(0)))
        pos = match.end()
    return out, data[pos:]


def _same_top_line_with_float_ulp(a: str, b: str, section: str | None) -> bool:
    if a == b:
        return True
    if section not in _TOP_FLOAT_SECTIONS:
        a_data, a_comment = _split_top_data_comment(a)
        b_data, b_comment = _split_top_data_comment(b)
        return a_comment == b_comment and a_data.split() == b_data.split()
    a_data, a_comment = _split_top_data_comment(a)
    b_data, b_comment = _split_top_data_comment(b)
    if a_comment != b_comment:
        return False
    a_tokens, a_trailing = _top_tokens_with_layout(a_data)
    b_tokens, b_trailing = _top_tokens_with_layout(b_data)
    if len(a_tokens) != len(b_tokens):
        return False
    for token_idx, ((a_sep, a_tok), (b_sep, b_tok)) in enumerate(zip(a_tokens, b_tokens)):
        if a_tok == b_tok:
            continue
        if not (_float_like(a_tok) and _float_like(b_tok)):
            return False
        a_float = _as_float(a_tok)
        b_float = _as_float(b_tok)
        if a_float is None or b_float is None:
            return False
        if (
            section == "dihedrals"
            and token_idx == 5
            and abs(abs(a_float) - 180.0) <= _TOP_FLOAT_ABS_TOL
            and abs(abs(b_float) - 180.0) <= _TOP_FLOAT_ABS_TOL
        ):
            continue
        limit = max(_TOP_FLOAT_ABS_TOL, _TOP_FLOAT_REL_TOL * max(abs(a_float), abs(b_float), 1.0))
        if abs(a_float - b_float) > limit:
            return False
    return True


def _top_matches_with_float_ulp(a_lines: list[str], b_lines: list[str]) -> bool:
    a_body = [line for line in _drop_top_header_metadata(a_lines) if line.strip()]
    b_body = [line for line in _drop_top_header_metadata(b_lines) if line.strip()]
    if len(a_body) != len(b_body):
        return False
    section: str | None = None
    for a_line, b_line in zip(a_body, b_body):
        a_section = _top_section(a_line)
        b_section = _top_section(b_line)
        if a_section or b_section:
            if a_section != b_section:
                return False
            section = a_section
            continue
        if not _same_top_line_with_float_ulp(a_line, b_line, section):
            return False
    return True


def _compare_file(a: Path, b: Path) -> dict:
    if not a.exists() and not b.exists():
        return {"match": True, "ignored": "both files absent"}
    if not a.exists() or not b.exists():
        return {"match": False, "reason": "missing file"}
    ta, tb = a.read_text().splitlines(), b.read_text().splitlines()
    if ta == tb:
        return {"match": True}
    if a.suffix == ".top" and _drop_top_header_metadata(ta) == _drop_top_header_metadata(tb):
        return {"match": True, "ignored": "topology header metadata before first section"}
    if a.suffix == ".xml" and _drop_xml_header_metadata(ta) == _drop_xml_header_metadata(tb):
        return {"match": True, "ignored": "OpenSMOG XML generated comment metadata before root element"}
    if a.suffix == ".top" and _top_matches_with_float_ulp(ta, tb):
        return {
            "match": True,
            "ignored": "topology header metadata, harmless whitespace layout, tiny floating-point print ULPs, and dihedral +/-180 endpoint print convention",
        }
    if max(len(ta), len(tb)) > 50000:
        for idx, (left, right) in enumerate(itertools.zip_longest(ta, tb, fillvalue="<missing>"), start=1):
            if left != right:
                return {
                    "match": False,
                    "diff": f"large-file first differing line {idx}\n--- {a}\n+++ {b}\n- {left}\n+ {right}",
                }
        return {"match": False, "diff": "large files differ"}
    diff_iter = difflib.unified_diff(ta, tb, fromfile=str(a), tofile=str(b), n=2)
    diff = "\n".join(itertools.islice(diff_iter, 120))
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
