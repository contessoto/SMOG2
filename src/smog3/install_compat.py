from __future__ import annotations

import argparse
import os
import stat
from pathlib import Path

TOOLS = ["adjustPDB", "tablegen", "extract", "scale-energies", "ions", "modifyXML", "editgro"]


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Install SMOG2-compatible shell wrappers for smog-check")
    p.add_argument("--bin-dir", default=".smog3-bin", help="Output directory for wrapper scripts")
    p.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]), help="Path to SMOG2 checkout")
    p.add_argument("--perl", default=os.environ.get("perl4smog", "perl"), help="Perl interpreter path")
    ns = p.parse_args(argv)

    repo = Path(ns.repo_root).resolve()
    bindir = Path(ns.bin_dir).resolve()
    bindir.mkdir(parents=True, exist_ok=True)

    header = (
        "#!/bin/bash\n"
        f"export PERLLIB={repo / 'src'}:$PERLLIB\n"
        f"export PERL5LIB={repo / 'src'}:$PERL5LIB\n"
        f"export SMOG_PATH={repo}\n"
        f"export perl4smog={ns.perl}\n"
    )

    smog2 = header + f"{ns.perl} {repo / 'src' / 'smogv2'} \"$@\"\n"
    _write_executable(bindir / "smog2", smog2)

    for tool in TOOLS:
        script = header + f"{ns.perl} {repo / 'src' / 'tools' / tool} \"$@\"\n"
        _write_executable(bindir / f"smog_{tool}", script)

    print(f"Wrote compatibility wrappers to {bindir}")
    print(f"Add to PATH: export PATH={bindir}:$PATH")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
