from __future__ import annotations

import os
import pathlib
import subprocess
import sys
from typing import Sequence


ROOT = pathlib.Path(__file__).resolve().parents[2]


def _run_perl_script(script: pathlib.Path, argv: Sequence[str]) -> int:
    perl = os.environ.get("perl4smog", "perl")
    env = os.environ.copy()
    src_dir = str(ROOT / "src")
    env["PERLLIB"] = f"{src_dir}:{env.get('PERLLIB', '')}".rstrip(":")
    env["PERL5LIB"] = f"{src_dir}:{env.get('PERL5LIB', '')}".rstrip(":")
    env.setdefault("SMOG_PATH", str(ROOT))

    proc = subprocess.run([perl, str(script), *argv], env=env)
    return int(proc.returncode)


def run_smog2(argv: Sequence[str]) -> int:
    return _run_perl_script(ROOT / "src" / "smogv2", argv)


def run_tool(tool_name: str, argv: Sequence[str]) -> int:
    return _run_perl_script(ROOT / "src" / "tools" / tool_name, argv)
