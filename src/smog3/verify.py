from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from .install_compat import main as install_compat_main


def _perl_deps_ok() -> bool:
    probe = subprocess.run(
        ["perl", "-MXML::Simple", "-MXML::Validator::Schema", "-MPDL", "-e", "print 'ok'"],
        capture_output=True,
        text=True,
    )
    return probe.returncode == 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run smog3 verification checks and optional full smog-check")
    p.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    p.add_argument("--require-full", action="store_true", help="Fail if full smog-check cannot be executed")
    ns = p.parse_args(argv)

    repo = Path(ns.repo_root).resolve()
    env = {**os.environ, "PYTHONPATH": str(repo / "src")}

    pytest = subprocess.run(["python", "-m", "pytest", "-q", "tests"], cwd=repo, env=env)
    if pytest.returncode != 0:
        return pytest.returncode

    if not _perl_deps_ok():
        print("Perl deps missing for full smog-check (XML::Simple/XML::Validator::Schema/PDL).")
        return 3 if ns.require_full else 0

    bindir = repo / ".smog3-bin"
    if bindir.exists():
        shutil.rmtree(bindir)
    install_compat_main(["--bin-dir", str(bindir), "--repo-root", str(repo)])

    env2 = os.environ.copy()
    env2["PATH"] = f"{bindir}:{env2.get('PATH','')}"
    env2["PERLLIB"] = f"{repo / 'src'}:{env2.get('PERLLIB','')}"
    env2["PERL5LIB"] = f"{repo / 'src'}:{env2.get('PERL5LIB','')}"
    rc = subprocess.run(["./smog-check"], cwd=repo / "SMOG-CHECK", env=env2)
    return int(rc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
