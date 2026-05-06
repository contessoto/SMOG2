from __future__ import annotations

import argparse
from pathlib import Path

from .gmx import parse_top_sections, write_top_sections


def _parse_ions_def(template: Path, ionnm: str):
    ionf = template / "ions.def"
    for ln in ionf.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith(";"):
            continue
        p = s.split()
        if p[0] == ionnm:
            return float(p[1]), float(p[2]), float(p[3]), float(p[4])
    raise SystemExit(f"Could not find {ionnm} in ions.def file")


def _parse_gro(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[1].strip())
    return lines[0], lines[2:2+n], lines[2+n]


def _fmt_gro_atom(resi, resn, name, idx, x=0.0, y=0.0, z=0.0):
    return f"{resi:5d}{resn:<5}{name:>5}{idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-f", default="smog.top")
    p.add_argument("-g", default="smog.gro")
    p.add_argument("-of", default="smog.ions.top")
    p.add_argument("-og", default="smog.ions.gro")
    p.add_argument("-ionn", type=int, default=0)
    p.add_argument("-ionnm", default=None)
    p.add_argument("-ionq", type=float, default=None)
    p.add_argument("-ionm", type=float, default=None)
    p.add_argument("-ionC12", type=float, default=None)
    p.add_argument("-ionC6", type=float, default=None)
    p.add_argument("-t", default=None)
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_ions -f in.top -g in.gro -ionnm NAME -ionn N [-ionq Q -ionm M -ionC12 C12 -ionC6 C6 | -t template]")
        return 1
    if ns.ionnm is None or ns.ionn < 1:
        raise SystemExit("Please indicate ion name with -ionnm and positive number with -ionn")

    if ns.t:
        m, q, c12, c6 = _parse_ions_def(Path(ns.t), ns.ionnm)
    else:
        if None in (ns.ionm, ns.ionq, ns.ionC12, ns.ionC6):
            raise SystemExit("Must provide -ionm -ionq -ionC12 -ionC6 when -t is not used")
        m, q, c12, c6 = ns.ionm, ns.ionq, ns.ionC12, ns.ionC6

    sections = parse_top_sections(Path(ns.f))

    # atomtypes append
    for sec in sections:
        if sec.name == "atomtypes":
            sec.lines.append(f"{ns.ionnm}\t{m}\t{q}\tA\t{c6}\t{c12}\n")

    # add moleculetype/atoms before system
    insert_idx = next((i for i, s in enumerate(sections) if s.name == "system"), len(sections))
    from .gmx import TopSection
    sections.insert(insert_idx, TopSection(name="atoms", lines=[f"1 {ns.ionnm} 1 {ns.ionnm} {ns.ionnm} 1\n"]))
    sections.insert(insert_idx, TopSection(name="moleculetype", lines=[f"{ns.ionnm} 1\n"]))

    for sec in sections:
        if sec.name == "molecules":
            sec.lines.append(f"{ns.ionnm} {ns.ionn}\n")

    write_top_sections(Path(ns.of), sections)

    head, atoms, box = _parse_gro(Path(ns.g))
    start_idx = len(atoms) + 1
    # add ions at origin shifted by 0.01
    for i in range(ns.ionn):
        atoms.append(_fmt_gro_atom(99999, ns.ionnm[:5], ns.ionnm[:5], start_idx + i, 0.01 * i, 0.0, 0.0))
    out = [head, str(len(atoms)), *atoms, box]
    Path(ns.og).write_text("\n".join(out) + "\n")
    return 0
