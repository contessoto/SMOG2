from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from .gmx import parse_ndx, parse_top_sections, write_top_sections


@dataclass
class ScaleConfig:
    group_d: str
    group_c1: str
    group_c2: str
    rescale_c: float = 1.0
    rescale_d: float = 1.0


def _split_comment(line: str) -> tuple[str, str]:
    if ";" in line:
        a, b = line.split(";", 1)
        return a.rstrip(), ";" + b
    return line.rstrip(), ""


def _is_comment_or_blank(line: str) -> bool:
    s = line.strip()
    return not s or s.startswith(";")


def _parse_defaults_combrule(sections) -> int:
    for sec in sections:
        if sec.name == "defaults":
            for ln in sec.lines:
                left, _ = _split_comment(ln)
                if left.strip():
                    parts = left.split()
                    if len(parts) >= 2:
                        return int(parts[1])
    return 1


def _scale_dihedral_line(line: str, group: set[int], factor: float) -> str:
    if _is_comment_or_blank(line):
        return line
    left, comment = _split_comment(line)
    parts = left.split()
    if len(parts) < 8:
        return line
    try:
        a, b, c, d = (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]))
        fn = int(parts[4])
    except ValueError:
        return line
    if fn == 1 and {a, b, c, d}.issubset(group):
        if factor == 0:
            return f";  {line.rstrip()}  ;  removed by smog_scale-energies\n"
        parts[6] = str(float(parts[6]) * factor)
        return "\t".join(parts) + (f"\t{comment} ; scaled by {factor}" if comment else f" ; scaled by {factor}") + "\n"
    return line


def _pair_matches(i: int, j: int, g1: set[int], g2: set[int]) -> bool:
    return (i in g1 and j in g2) or (i in g2 and j in g1)


def _scale_pair_line(line: str, g1: set[int], g2: set[int], factor: float, combrule: int) -> str:
    if _is_comment_or_blank(line):
        return line
    left, comment = _split_comment(line)
    parts = left.split()
    if len(parts) < 5:
        return line
    try:
        i, j = int(parts[0]), int(parts[1])
    except ValueError:
        return line
    if _pair_matches(i, j, g1, g2):
        if factor == 0:
            return f";  {line.rstrip()}  ;  removed by smog_scale-energies\n"
        if combrule == 1:
            parts[3] = str(float(parts[3]) * factor)
        parts[4] = str(float(parts[4]) * factor)
        return "\t".join(parts) + (f"\t{comment} ; scaled by {factor}" if comment else f" ; scaled by {factor}") + "\n"
    return line


def _scale_exclusion_line(line: str, g1: set[int], g2: set[int], factor: float) -> str:
    if factor != 0 or _is_comment_or_blank(line):
        return line
    left, _ = _split_comment(line)
    parts = left.split()
    if len(parts) < 2:
        return line
    try:
        i, j = int(parts[0]), int(parts[1])
    except ValueError:
        return line
    if _pair_matches(i, j, g1, g2):
        return f";  {line.rstrip()}  ; removed by smog_scale-energies\n"
    return line


def scale_topology(top_in: str | Path, ndx_in: str | Path, top_out: str | Path, cfg: ScaleConfig) -> None:
    groups = parse_ndx(ndx_in)
    dg = set(groups[cfg.group_d])
    c1 = set(groups[cfg.group_c1])
    c2 = set(groups[cfg.group_c2])

    sections = parse_top_sections(top_in)
    combrule = _parse_defaults_combrule(sections)

    for sec in sections:
        if sec.name == "dihedrals" and cfg.rescale_d != 1.0:
            sec.lines = [_scale_dihedral_line(ln, dg, cfg.rescale_d) for ln in sec.lines]
        elif sec.name == "pairs" and cfg.rescale_c != 1.0:
            sec.lines = [_scale_pair_line(ln, c1, c2, cfg.rescale_c, combrule) for ln in sec.lines]
        elif sec.name == "exclusions" and cfg.rescale_c == 0:
            sec.lines = [_scale_exclusion_line(ln, c1, c2, cfg.rescale_c) for ln in sec.lines]

    write_top_sections(top_out, sections)


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-f", default="smog.top")
    p.add_argument("-n", default="smog.ndx")
    p.add_argument("-of", default="smog.rescaled.top")
    p.add_argument("-rc", type=float, default=1.0)
    p.add_argument("-rd", type=float, default=1.0)
    p.add_argument("-grpD", default=None)
    p.add_argument("-grpC1", default=None)
    p.add_argument("-grpC2", default=None)
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_scale-energies -f smog.top -n smog.ndx -of out.top -rc 1.0 -rd 1.0 [-grpD name] [-grpC1 name -grpC2 name]")
        return 1

    if ns.f == ns.of:
        raise SystemExit(f"Input and output top files can not have the same name: {ns.of}")
    if ns.rc < 0 or ns.rd < 0:
        raise SystemExit("negative may not be provided with -rc/-rd")
    if ns.rc == 1.0 and ns.rd == 1.0:
        raise SystemExit("-rc and -rd given values of 1. No .top file generated.")

    groups = parse_ndx(ns.n)
    names = list(groups.keys())
    if not names:
        raise SystemExit("no atom groups given in ndx file.")

    grpD = ns.grpD or names[0]
    grpC1 = ns.grpC1 or names[0]
    grpC2 = ns.grpC2 or names[min(1, len(names)-1)]

    cfg = ScaleConfig(group_d=grpD, group_c1=grpC1, group_c2=grpC2, rescale_c=ns.rc, rescale_d=ns.rd)
    scale_topology(ns.f, ns.n, ns.of, cfg)
    print("\n\tSUCCESS: Interactions rescaled.\n")
    return 0
