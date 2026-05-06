from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class FileDigest:
    path: str
    sha256: str
    size: int


def digest_files(root: Path, patterns: Iterable[str]) -> list[FileDigest]:
    files: dict[str, FileDigest] = {}
    for pat in patterns:
        for p in root.glob(pat):
            if p.is_file():
                rel = str(p.relative_to(root))
                data = p.read_bytes()
                files[rel] = FileDigest(path=rel, sha256=hashlib.sha256(data).hexdigest(), size=len(data))
    return [files[k] for k in sorted(files)]


def compare_digest_sets(a: list[FileDigest], b: list[FileDigest]) -> dict:
    ad = {x.path: x for x in a}
    bd = {x.path: x for x in b}
    only_a = sorted(set(ad) - set(bd))
    only_b = sorted(set(bd) - set(ad))
    mismatched = []
    for k in sorted(set(ad) & set(bd)):
        if ad[k].sha256 != bd[k].sha256:
            mismatched.append(
                {
                    "path": k,
                    "a_sha256": ad[k].sha256,
                    "b_sha256": bd[k].sha256,
                    "a_size": ad[k].size,
                    "b_size": bd[k].size,
                }
            )
    return {
        "only_in_a": only_a,
        "only_in_b": only_b,
        "mismatched": mismatched,
        "match": not only_a and not only_b and not mismatched,
    }


def write_report(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
