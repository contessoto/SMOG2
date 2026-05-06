from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .install_compat import main as install_compat_main
from .parity_report import compare_digest_sets, digest_files, write_report


ARTIFACT_PATTERNS = ["FAILED/**/*", "FAILED.tools/**/*", "*.log", "src/*.log"]


def _run_smogcheck(repo_root: Path, bin_dir: Path, start: int, end: int, workdir: Path) -> int:
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env.get('PATH','')}"
    env["PERLLIB"] = f"{repo_root / 'src'}:{env.get('PERLLIB','')}"
    env["PERL5LIB"] = f"{repo_root / 'src'}:{env.get('PERL5LIB','')}"
    proc = subprocess.run(["./smog-check", str(start), str(end)], cwd=workdir, env=env)
    return int(proc.returncode)


def _prepare_workdir(repo_root: Path, label: str, temp_root: Path) -> Path:
    dest = temp_root / label
    shutil.copytree(repo_root / "SMOG-CHECK", dest)
    return dest


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run smog-check subset with smog3 wrappers and optional baseline comparison")
    p.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    p.add_argument("--start", type=int, default=1)
    p.add_argument("--end", type=int, default=1)
    p.add_argument("--perl", default=os.environ.get("perl4smog", "perl"))
    p.add_argument("--baseline-smog2", default=None, help="Path to baseline smog2 executable for two-run parity report")
    p.add_argument("--report-json", default="smog3_parity_report.json")
    ns = p.parse_args(argv)

    repo_root = Path(ns.repo_root).resolve()
    with tempfile.TemporaryDirectory(prefix="smog3-check-") as td:
        troot = Path(td)
        candidate_bin = troot / "candidate-bin"
        install_compat_main(["--bin-dir", str(candidate_bin), "--repo-root", str(repo_root), "--perl", ns.perl])

        candidate_work = _prepare_workdir(repo_root, "candidate", troot)
        rc = _run_smogcheck(repo_root, candidate_bin, ns.start, ns.end, candidate_work)
        if rc != 0:
            print(f"smog-check failed for tests {ns.start}..{ns.end} with exit code {rc}")
            return rc
        print(f"smog-check passed for tests {ns.start}..{ns.end} using smog3 compatibility wrappers")

        if ns.baseline_smog2:
            baseline_bin = troot / "baseline-bin"
            baseline_bin.mkdir(parents=True, exist_ok=True)
            shutil.copy2(Path(ns.baseline_smog2), baseline_bin / "smog2")
            (baseline_bin / "smog2").chmod(0o755)

            baseline_work = _prepare_workdir(repo_root, "baseline", troot)
            brc = _run_smogcheck(repo_root, baseline_bin, ns.start, ns.end, baseline_work)
            if brc != 0:
                print(f"Baseline smog-check failed for tests {ns.start}..{ns.end} with exit code {brc}")
                return brc

            cand = digest_files(candidate_work, ARTIFACT_PATTERNS)
            base = digest_files(baseline_work, ARTIFACT_PATTERNS)
            report = compare_digest_sets(base, cand)
            report["start"] = ns.start
            report["end"] = ns.end
            report["baseline_smog2"] = str(ns.baseline_smog2)
            write_report(Path(ns.report_json).resolve(), report)
            print(f"Wrote parity report to {Path(ns.report_json).resolve()}")
            if not report["match"]:
                return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
