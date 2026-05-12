from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class MapRules:
    """The subset of SMOG map-file behavior needed by native adjustPDB.

    SMOG tutorial map files use two styles that matter for model-generation
    parity: direct atom mappings, and residue declarations with terminal
    markers such as ``%first``/``%last``.  The latter lets custom templates
    rename ASN -> ASNN or GLY -> GLYC before smog2/smog3 assigns parameters.
    """

    atom_mapping: dict[tuple[str, str], tuple[str, str]] = field(default_factory=dict)
    global_atom_renames: dict[str, str] = field(default_factory=dict)
    first_residue_names: dict[str, str] = field(default_factory=dict)
    last_residue_names: dict[str, str] = field(default_factory=dict)


def _base_residue_name(name: str) -> str:
    if len(name) <= 3:
        return name
    if name.endswith(("N", "C", "T")):
        return name[:-1]
    return name[:3]


def _parse_map(path: Path) -> MapRules:
    rules = MapRules()
    if not path.exists():
        return rules
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith(";") or s.startswith("#"):
            continue
        parts = s.split()
        keyword = parts[0].lower()
        if keyword == "rename" and len(parts) >= 3:
            rules.global_atom_renames[parts[1]] = parts[2]
            continue
        if keyword == "residue" and len(parts) >= 2:
            name = parts[1]
            base = _base_residue_name(name)
            if "%first" in parts:
                rules.first_residue_names[base] = name
            if "%last" in parts:
                rules.last_residue_names[base] = name
            continue
        if len(parts) >= 4:
            rules.atom_mapping[(parts[0], parts[1])] = (parts[2], parts[3])
    return rules


def _read_insert_ter_answer() -> str:
    answer = sys.stdin.readline()
    return answer.strip().lower()


def _should_insert_ter(previous: tuple[str, int] | None, current: tuple[str, int], insert_all: bool) -> tuple[bool, bool]:
    if previous is None:
        return False, insert_all
    prev_chain, prev_resi = previous
    chain, resi = current
    nonsequential = chain != prev_chain or resi != prev_resi + 1
    if not nonsequential:
        return False, insert_all
    if insert_all:
        return True, insert_all
    answer = _read_insert_ter_answer()
    if answer.startswith("a"):
        return True, True
    return answer.startswith("y"), insert_all


def _fmt_atom(serial: int, name: str, altloc: str, resn: str, chain: str, resi: int, x: float, y: float, z: float, elem: str, rec: str="ATOM") -> str:
    return f"{rec:<6}{serial:>5d} {name:>4}{altloc:1}{resn:<4}{chain:1}{resi:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {elem:>2}"


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i", required=False)
    p.add_argument("-o", default="adjusted.pdb")
    p.add_argument("-map", default=None)
    p.add_argument("-default", action="store_true")
    p.add_argument("-legacy", action="store_true")
    p.add_argument("-altloc", action="store_true")
    p.add_argument("-large", action="store_true")
    p.add_argument("-removewater", action="store_true")
    p.add_argument("-PDBnums", action="store_true")
    p.add_argument("-removeH", action="store_true")
    p.add_argument("-sort", action="store_true")
    p.add_argument("-insertTER", action="store_true")
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra or not ns.i:
        print("usage: smog_adjustPDB -i input.pdb [-o adjusted.pdb] [-map file] [-altloc] [-removeH] [-removewater] [-PDBnums] [-sort] [-insertTER] [-large]")
        return 1

    rules = _parse_map(Path(ns.map)) if ns.map else MapRules()

    recs = []
    serial_out = 0
    cur_res = None
    res_out = 0
    lines = Path(ns.i).read_text().splitlines()
    previous_residue: tuple[str, int] | None = None
    insert_all_ter = False
    pending_ter = False
    for ln in lines:
        if ln.startswith("TER"):
            pending_ter = True
            previous_residue = None
            continue
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        rec = ln[:6].strip()
        serial = int(ln[6:11])
        name = ln[12:16].strip()
        altloc = ln[16:17]
        resn = ln[17:20].strip()
        chain = ln[21:22]
        resi = int(ln[22:26])
        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
        elem = (ln[76:78].strip() or name[:1]).upper()

        if ns.removewater and resn in {"HOH", "WAT", "SOL"}:
            continue
        if ns.removeH and (elem == "H" or name.startswith("H")):
            continue
        if not ns.altloc and altloc not in {" ", "A"}:
            continue

        original_residue = (chain, resi)
        ter_before = pending_ter
        reset_before = pending_ter
        pending_ter = False
        if ns.insertTER and not ter_before and previous_residue != original_residue:
            ter_before, insert_all_ter = _should_insert_ter(previous_residue, original_residue, insert_all_ter)
            reset_before = False
        previous_residue = original_residue

        name = rules.global_atom_renames.get(name, name)
        key = (resn, name)
        if key in rules.atom_mapping:
            resn, name = rules.atom_mapping[key]

        if ns.sort:
            sort_key = (chain, resi, name, serial)
        else:
            sort_key = (len(recs),)

        recs.append([sort_key, rec, serial, name, altloc if ns.altloc else " ", resn, chain, resi, x, y, z, elem, ter_before, reset_before])

    recs.sort(key=lambda t: t[0])

    residue_order: list[tuple[str, int, int]] = []
    segment = 0
    seen_residues: set[tuple[str, int, int]] = set()
    for row in recs:
        if row[12]:
            segment += 1
        key = (str(row[6]), int(row[7]), segment)
        if key not in seen_residues:
            residue_order.append(key)
            seen_residues.add(key)
    first_residues = set()
    last_residues = set()
    by_segment: dict[int, list[tuple[str, int, int]]] = {}
    for key in residue_order:
        by_segment.setdefault(key[2], []).append(key)
    for items in by_segment.values():
        if items:
            first_residues.add(items[0])
            last_residues.add(items[-1])

    out_lines = []
    if ns.large:
        out_lines.append("LARGE FORMAT")
    out_lines.append("REMARK adjusted by smog_adjustPDB (python-native)")

    segment = 0
    for _, rec, serial, name, altloc, resn, chain, resi, x, y, z, elem, ter_before, reset_before in recs:
        if ter_before:
            out_lines.append("TER")
            segment += 1
            cur_res = None
            if reset_before:
                serial_out = 0
                res_out = 0
        residue_key = (chain, resi, segment)
        base = _base_residue_name(resn)
        if residue_key in first_residues and base in rules.first_residue_names:
            resn = rules.first_residue_names[base]
        if residue_key in last_residues and base in rules.last_residue_names:
            resn = rules.last_residue_names[base]
        serial_out = serial if ns.PDBnums else serial_out + 1
        if ns.PDBnums:
            res_out = resi
        else:
            if cur_res != (chain, resi, segment):
                res_out += 1
                cur_res = (chain, resi, segment)
        out_lines.append(_fmt_atom(serial_out, name, altloc, resn, chain, res_out, x, y, z, elem, rec))

    Path(ns.o).write_text("\n".join(out_lines) + "\n")
    return 0
