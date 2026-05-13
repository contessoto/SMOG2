"""Python-native subset of SMOG ``smog_extract``.

This runtime helper extracts atom subsets from GRO, NDX, topology, contacts, and
OpenSMOG XML files while remapping atom indices consistently.  It is used by
large-fragment tutorial workflows and does not call Perl or original SMOG2.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

from .gmx import parse_ndx, parse_top_sections, write_top_sections


def _parse_gro(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0]
    natoms = int(lines[1].strip())
    atoms = lines[2 : 2 + natoms]
    box = lines[2 + natoms]
    return header, atoms, box


def _write_gro(path: Path, header: str, atoms: list[str], box: str):
    """Write an extracted GRO file with the retained atom lines and box."""

    out = [header, str(len(atoms)), *atoms, box]
    path.write_text("\n".join(out) + "\n")


def _select_indices(args, groups, atoms):
    if args.nondx:
        if args.distfrom is None or args.distval is None:
            raise SystemExit("-nondx currently requires -distfrom and -distval")
        ref = args.distfrom
        if ref < 1 or ref > len(atoms):
            raise SystemExit("-distfrom index out of range")
        xc = float(atoms[ref - 1][20:28])
        yc = float(atoms[ref - 1][28:36])
        zc = float(atoms[ref - 1][36:44])
        d2 = args.distval * args.distval
        keep = []
        for i, ln in enumerate(atoms, start=1):
            x = float(ln[20:28]); y = float(ln[28:36]); z = float(ln[36:44])
            if (x - xc) ** 2 + (y - yc) ** 2 + (z - zc) ** 2 < d2:
                keep.append(i)
        return keep

    names = list(groups.keys())
    if not names:
        raise SystemExit("no atom groups found in ndx file")
    group = args.group if args.group else names[0]
    if group not in groups:
        raise SystemExit(f"group {group} not found")
    keep = groups[group]
    return keep if args.ndxorder else sorted(keep)


def _remap_interaction_line(line: str, remap: dict[int, int], nidx: int):
    s = line.strip()
    if not s or s.startswith(";"):
        return line
    parts = line.split()
    if len(parts) < nidx:
        return None
    try:
        ids = [int(parts[i]) for i in range(nidx)]
    except ValueError:
        return line
    if not all(i in remap for i in ids):
        return None
    for i in range(nidx):
        parts[i] = str(remap[ids[i]])
    return "\t".join(parts) + "\n"


def _write_filtered_opensmog_xml(src: Path, dst: Path, remap: dict[int, int]) -> None:
    """Filter OpenSMOG interactions to the extracted atom set.

    The legacy ``smog_extract`` tool preserves the OpenSMOG force blocks but
    drops interactions whose atom indices were removed, remapping retained
    indices to the extracted topology.  A line-preserving implementation keeps
    SMOG3 output comparable to SMOG2 while avoiding a full XML reformatter.
    """

    atom_attr = re.compile(r'\b([ijkl])="(\d+)"')
    out: list[str] = []
    for line in src.read_text().splitlines():
        if "<interaction " not in line:
            out.append(line)
            continue
        ids = {name: int(value) for name, value in atom_attr.findall(line)}
        if not ids or not all(value in remap for value in ids.values()):
            continue

        def repl(match: re.Match[str]) -> str:
            return f'{match.group(1)}="{remap[int(match.group(2))]}"'

        out.append(atom_attr.sub(repl, line))
    dst.write_text("\n".join(out) + "\n")


def main(argv: list[str]) -> int:
    """Run native extraction for selected NDX groups or distance selections."""

    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-f", default="smog.top")
    p.add_argument("-g", default="smog.gro")
    p.add_argument("-n", default="smog.ndx")
    p.add_argument("-of", default="extracted.top")
    p.add_argument("-og", default="extracted.gro")
    p.add_argument("-om", default="atomindex.map")
    p.add_argument("-orm", default="restrained.map")
    p.add_argument("-group", default=None)
    p.add_argument("-ndxorder", action="store_true")
    p.add_argument("-nondx", action="store_true")
    p.add_argument("-distval", type=float, default=None)
    p.add_argument("-distfrom", type=int, default=None)
    p.add_argument("-help", "-?", action="store_true")
    p.add_argument("-restraints", default="null")
    p.add_argument("-OpenSMOG", default=None)
    p.add_argument("-OpenSMOGout", default="extracted.xml")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_extract -f in.top -g in.gro -n in.ndx -of out.top -og out.gro -om atomindex.map")
        return 1
    opensmog_in = Path(ns.OpenSMOG) if ns.OpenSMOG else None

    top_in = Path(ns.f); gro_in = Path(ns.g); ndx_in = Path(ns.n)
    out_top = Path(ns.of); out_gro = Path(ns.og); out_map = Path(ns.om)

    header, atoms, box = _parse_gro(gro_in)
    groups = {} if ns.nondx else parse_ndx(ndx_in)
    keep = _select_indices(ns, groups, atoms)

    remap = {old: i + 1 for i, old in enumerate(keep)}

    selected_atoms = [atoms[i - 1] for i in keep]
    _write_gro(out_gro, header + ". Note: This file was processed by smog_extract (version 2.4.5)", selected_atoms, box)

    with out_map.open("w") as fh:
        fh.write("This file contains the corresponding atom indices for the new (left column) and old (right column) systems.\n")
        for new, old in sorted((v, k) for k, v in remap.items()):
            fh.write(f"{new} {old}\n")

    sections = parse_top_sections(top_in)
    for sec in sections:
        if sec.name == "atoms":
            new_lines = []
            for ln in sec.lines:
                if ln.strip().startswith(";") or not ln.strip():
                    new_lines.append(ln); continue
                parts = ln.split()
                if not parts: continue
                try:
                    idx = int(parts[0])
                except ValueError:
                    new_lines.append(ln); continue
                if idx in remap:
                    parts[0] = str(remap[idx]); parts[5] = str(remap[idx])
                    new_lines.append("\t".join(parts) + "\n")
            sec.lines = new_lines
        elif sec.name in {"bonds", "pairs", "exclusions"}:
            sec.lines = [x for x in (_remap_interaction_line(ln, remap, 2) for ln in sec.lines) if x is not None]
        elif sec.name in {"angles"}:
            sec.lines = [x for x in (_remap_interaction_line(ln, remap, 3) for ln in sec.lines) if x is not None]
        elif sec.name in {"dihedrals"}:
            sec.lines = [x for x in (_remap_interaction_line(ln, remap, 4) for ln in sec.lines) if x is not None]
        elif sec.name == "defaults":
            normalized = []
            for ln in sec.lines:
                s = ln.strip()
                if not s:
                    continue
                if s.startswith(";"):
                    normalized.append(ln)
                    continue
                normalized.append(" ".join(s.split()) + "\n")
            sec.lines = normalized

    write_top_sections(out_top, sections)

    if ns.restraints != "null":
        Path(ns.orm).write_text("# restrained atoms generated by smog_extract\n" + "\n".join(str(v) for v in sorted(remap.values())) + "\n")

    if opensmog_in is not None:
        text = opensmog_in.read_text()
        if "<OpenSMOGforces" in text and "<interaction " in text:
            _write_filtered_opensmog_xml(opensmog_in, Path(ns.OpenSMOGout), remap)
        else:
            Path(ns.OpenSMOGout).write_text(text)

    print("\n\tSUCCESS: Extraction complete.\n")
    return 0


if __name__ == "__main__":
    import sys

    raise SystemExit(main(sys.argv[1:]))
