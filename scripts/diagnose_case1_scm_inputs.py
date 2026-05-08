#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "parity_runs" / "case1" / "baseline"
CAND = ROOT / "parity_runs" / "case1" / "candidate"


def _gro_atom_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        return []
    try:
        count = int(lines[1].strip())
    except ValueError:
        return []
    return lines[2:2 + count]


def _gro_summary(path: Path) -> dict:
    atom_lines = _gro_atom_lines(path)
    coords = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) >= 6:
            coords.append(tuple(float(x) for x in parts[-3:]))
    summary = {
        "exists": path.exists(),
        "atom_count": len(atom_lines),
        "first20": atom_lines[:20],
        "last20": atom_lines[-20:],
    }
    if coords:
        summary["coord_min"] = [min(axis) for axis in zip(*coords)]
        summary["coord_max"] = [max(axis) for axis in zip(*coords)]
    return summary


def _ndx_summary(path: Path) -> dict:
    groups = []
    current = None
    for raw in path.read_text().splitlines() if path.exists() else []:
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = {"name": s.strip("[] "), "atoms": 0}
            groups.append(current)
        elif current and s:
            current["atoms"] += len(s.split())
    return {"exists": path.exists(), "group_count": len(groups), "groups": groups}


def _top_summary(path: Path) -> dict:
    counts = {"atoms": 0, "bonds": 0, "angles": 0, "dihedrals": 0}
    current = None
    for raw in path.read_text().splitlines() if path.exists() else []:
        s = raw.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s.strip("[] ")
            continue
        if current in counts and s and not s.startswith(";"):
            counts[current] += 1
    return {"exists": path.exists(), **counts}


def _contacts(path: Path) -> set[tuple[str, str, str, str]]:
    out = set()
    for raw in path.read_text().splitlines() if path.exists() else []:
        parts = raw.split()
        if len(parts) >= 4:
            out.add(tuple(parts[:4]))
    return out


def _chain_pair_counts(contacts: set[tuple[str, str, str, str]]) -> dict[str, int]:
    counts = Counter(f"{a}-{c}" for a, _b, c, _d in contacts)
    return dict(sorted(counts.items()))


def main() -> int:
    baseline_contacts = _contacts(BASE / "model.contacts")
    candidate_contacts = _contacts(CAND / "model.contacts")
    only_baseline = sorted(baseline_contacts - candidate_contacts)
    only_candidate = sorted(candidate_contacts - baseline_contacts)
    report = {
        "baseline": {
            "model.gro": _gro_summary(BASE / "model.gro"),
            "model.gro4SCM.gro": _gro_summary(BASE / "model.gro4SCM.gro"),
            "model.ndx": _ndx_summary(BASE / "model.ndx"),
            "model.top4SCM.top": _top_summary(BASE / "model.top4SCM.top"),
            "scm_command": "java -jar SCM.jar -g model.gro4SCM.gro -freecoor -t model.top4SCM.top -o model.contacts -ch model.ndx -m shadow -c 6 -s 1 -br 0.5 -bif AA-whitford09.bif -pd 3 --smog2output --showProgress",
        },
        "candidate": {
            "model.gro": _gro_summary(CAND / "model.gro"),
            "model.gro4SCM.gro": _gro_summary(CAND / "model.gro4SCM.gro"),
            "model.ndx": _ndx_summary(CAND / "model.ndx"),
            "model.top4SCM.top": _top_summary(CAND / "model.top4SCM.top"),
            "scm_command": "java -jar SCM.jar -g model.gro4SCM.gro -freecoor -t model.top4SCM.top -o model.contacts -ch model.ndx -m shadow -c 6 -s 1 -br 0.5 -bif AA-whitford09.bif -pd 3 --smog2output --showProgress",
        },
        "contacts": {
            "baseline_count": len(baseline_contacts),
            "candidate_count": len(candidate_contacts),
            "only_baseline_count": len(only_baseline),
            "only_candidate_count": len(only_candidate),
            "only_baseline_first50": [" ".join(x) for x in only_baseline[:50]],
            "only_candidate_first50": [" ".join(x) for x in only_candidate[:50]],
            "baseline_chain_pairs": _chain_pair_counts(baseline_contacts),
            "candidate_chain_pairs": _chain_pair_counts(candidate_contacts),
        },
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
