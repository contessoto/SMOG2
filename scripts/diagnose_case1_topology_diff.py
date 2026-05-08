#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "parity_runs" / "case1" / "baseline" / "model.top"
CAND = ROOT / "parity_runs" / "case1" / "candidate" / "model.top"

KEY_COLUMNS = {
    "bonds": 2,
    "angles": 3,
    "dihedrals": 5,
    "pairs": 2,
    "exclusions": 2,
}


def _sections(path: Path) -> dict[str, list[str]]:
    sections: dict[str, list[str]] = {"comments": []}
    current = "comments"
    for raw in path.read_text().splitlines() if path.exists() else []:
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            sections.setdefault(current, [])
            continue
        if not s:
            continue
        sections.setdefault(current, []).append(raw.rstrip())
    return sections


def _content(lines: list[str]) -> list[str]:
    return [line.split(";", 1)[0].strip() for line in lines if line.split(";", 1)[0].strip()]


def _first_diff(left: list[str], right: list[str]) -> dict[str, str | int] | None:
    for idx, (a, b) in enumerate(zip(left, right), start=1):
        if a != b:
            return {"line": idx, "baseline": a, "candidate": b}
    if len(left) != len(right):
        return {"line": min(len(left), len(right)) + 1, "baseline": "<EOF>" if len(left) < len(right) else left[-1], "candidate": "<EOF>" if len(right) < len(left) else right[-1]}
    return None


def _classify(section: str, base: list[str], cand: list[str]) -> str:
    if base == cand:
        return "exact"
    base_content = _content(base)
    cand_content = _content(cand)
    if base_content == cand_content:
        return "metadata/comment-or-whitespace only"
    if Counter(base_content) == Counter(cand_content):
        return "ordering or formatting"
    key_cols = KEY_COLUMNS.get(section)
    if key_cols is not None:
        base_keys = Counter(tuple(line.split()[:key_cols]) for line in base_content)
        cand_keys = Counter(tuple(line.split()[:key_cols]) for line in cand_content)
        if base_keys == cand_keys:
            return "numerical parameter value or per-line formatting"
        base_canonical = Counter(tuple(sorted(line.split()[:key_cols])) for line in base_content)
        cand_canonical = Counter(tuple(sorted(line.split()[:key_cols])) for line in cand_content)
        if base_canonical == cand_canonical:
            return "ordering/orientation"
    return "missing/wrong scientific value"


def main() -> int:
    baseline = _sections(BASE)
    candidate = _sections(CAND)
    names = [name for name in baseline.keys() if name in candidate]
    for name in candidate:
        if name not in baseline:
            names.append(name)
    report = {}
    for name in names:
        base = baseline.get(name, [])
        cand = candidate.get(name, [])
        base_content = _content(base)
        cand_content = _content(cand)
        report[name] = {
            "classification": _classify(name, base, cand),
            "baseline_lines": len(base_content),
            "candidate_lines": len(cand_content),
            "ordered_match": base == cand,
            "content_match": base_content == cand_content,
            "content_multiset_match": Counter(base_content) == Counter(cand_content),
            "first_diff": _first_diff(base, cand),
        }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
