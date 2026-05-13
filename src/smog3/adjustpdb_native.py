"""Python-native subset of SMOG ``smog_adjustPDB`` preprocessing.

The public tutorial workflows use ``smog_adjustPDB`` to remove waters, strip
hydrogens, insert TER records, normalize nucleic atom names, apply map files,
and rename terminal/custom residues before model generation.  This runtime
module implements those preprocessing behaviors directly in Python and never
calls the legacy Perl tool.
"""

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


_NUCLEIC_RESIDUES = {
    "A", "C", "G", "U", "I",
    "DA", "DC", "DG", "DT", "DI",
    "ADE", "CYT", "GUA", "URA", "THY",
}

_AMINO_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

_DEFAULT_RESIDUE_RENAMES = {
    "1MG": "G",
    "2MG": "G",
    "OMG": "G",
    "MEQ": "ALA",
    "2MA": "A",
    "MA6": "A",
    "4OC": "C",
    "5MC": "C",
    "OMC": "C",
    "5MU": "U",
    "D2T": "ASP",
    "H2U": "U",
    "OMU": "U",
    "UR3": "U",
}


def _smog2_default_residue_name(resn: str) -> str:
    return _DEFAULT_RESIDUE_RENAMES.get(resn.upper(), resn)


def _smog2_default_atom_name(resn: str, atom_name: str) -> str:
    """Apply the built-in SMOG2 nucleic-acid atom-name compatibility map.

    SMOG2's ``smog_adjustPDB`` converts common PDB v3 nucleic names such as
    ``OP1``/``O5'`` to the template-era names ``O1P``/``O5*`` even when no
    explicit ``-map`` file is supplied.  Large-system and tutorial workflows
    depend on those names because the AA templates encode nucleic geometry with
    the star/phosphate aliases.
    """

    base = _base_residue_name(resn).upper()
    if base not in _NUCLEIC_RESIDUES:
        return atom_name
    phosphate = {"OP1": "O1P", "OP2": "O2P", "OP3": "O3P"}
    if atom_name in phosphate:
        return phosphate[atom_name]
    return atom_name.replace("'", "*")


def _smog2_default_terminal_residue_name(
    resn: str,
    atom_names: set[str],
    *,
    first_residue: bool,
    last_residue: bool,
) -> str:
    """Apply SMOG2's built-in terminal residue aliases for adjusted PDBs."""

    base = _base_residue_name(resn).upper()
    if last_residue and base in _AMINO_RESIDUES and "OXT" in atom_names:
        return f"{base}T"
    if first_residue and base in _NUCLEIC_RESIDUES and "P" not in atom_names:
        return f"{base}0P"
    if base in _AMINO_RESIDUES and atom_names == {"N", "CA", "C", "O", "CB"}:
        return "ALA"
    return resn


def _read_insert_ter_answer() -> str:
    answer = sys.stdin.readline()
    return answer.strip().lower()


def _should_insert_ter(previous: tuple[str, int, str] | None, current: tuple[str, int, str], insert_all: bool) -> tuple[bool, bool]:
    if previous is None:
        return False, insert_all
    prev_chain, prev_resi, _prev_icode = previous
    chain, resi, _icode = current
    nonsequential = chain != prev_chain or (resi != prev_resi and resi != prev_resi + 1)
    if not nonsequential:
        return False, insert_all
    if insert_all:
        return True, insert_all
    answer = _read_insert_ter_answer()
    if answer.startswith("a"):
        return True, True
    return answer.startswith("y"), insert_all


def _fmt_atom(serial: int, name: str, altloc: str, resn: str, chain: str, resi: int, icode: str, x: float, y: float, z: float, elem: str, rec: str="ATOM") -> str:
    return f"{rec:<6}{serial:>5d} {name:>4}{altloc:1}{resn:<4}{chain:1}{resi:>4d}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {elem:>2}"


def main(argv: list[str]) -> int:
    """Run native PDB adjustment and write the adjusted PDB to disk.

    Supported options cover the tutorial and SMOG-CHECK workflows: water/H
    removal, map-file residue/atom renaming, large-format numbering, insertion
    of TER records, and SMOG2-compatible default residue normalization.
    """

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
    previous_residue: tuple[str, int, str] | None = None
    insert_all_ter = False
    pending_ter = False
    pending_ter_after_skipped = False
    skipped_since_output = False
    for ln in lines:
        if ln.startswith("TER"):
            pending_ter = True
            pending_ter_after_skipped = skipped_since_output
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
        icode = ln[26:27]
        x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
        elem = (ln[76:78].strip() or name[:1]).upper()

        if ns.removewater and resn in {"HOH", "WAT", "SOL"}:
            skipped_since_output = True
            continue
        if ns.removeH and (elem == "H" or name.startswith("H")):
            skipped_since_output = True
            continue
        if not ns.altloc and altloc not in {" ", "A"}:
            skipped_since_output = True
            continue

        original_residue = (chain, resi, icode)
        ter_before = 1 if pending_ter else 0
        if pending_ter and pending_ter_after_skipped and ns.insertTER:
            ter_before = 2
        reset_before = pending_ter
        pending_ter = False
        pending_ter_after_skipped = False
        if ns.insertTER and not ter_before and previous_residue != original_residue:
            ter_before, insert_all_ter = _should_insert_ter(previous_residue, original_residue, insert_all_ter)
            ter_before = 1 if ter_before else 0
            reset_before = False
        previous_residue = original_residue
        skipped_since_output = False

        if ns.map is None or ns.default:
            resn = _smog2_default_residue_name(resn)
            name = _smog2_default_atom_name(resn, name)
        name = rules.global_atom_renames.get(name, name)
        key = (resn, name)
        if key in rules.atom_mapping:
            resn, name = rules.atom_mapping[key]

        if ns.sort:
            sort_key = (chain, resi, name, serial)
        else:
            sort_key = (len(recs),)

        recs.append([sort_key, rec, serial, name, altloc if ns.altloc else " ", resn, chain, resi, icode, x, y, z, elem, ter_before, reset_before])

    recs.sort(key=lambda t: t[0])

    residue_order: list[tuple[str, int, str, int]] = []
    segment = 0
    seen_residues: set[tuple[str, int, str, int]] = set()
    for row in recs:
        if row[13]:
            segment += int(row[13])
        key = (str(row[6]), int(row[7]), str(row[8]), segment)
        if key not in seen_residues:
            residue_order.append(key)
            seen_residues.add(key)
    first_residues = set()
    last_residues = set()
    by_segment: dict[int, list[tuple[str, int, str, int]]] = {}
    for key in residue_order:
        by_segment.setdefault(key[3], []).append(key)
    for items in by_segment.values():
        if items:
            first_residues.add(items[0])
            last_residues.add(items[-1])
    atom_names_by_residue: dict[tuple[str, int, str, int], set[str]] = {}
    segment = 0
    for row in recs:
        if row[13]:
            segment += int(row[13])
        residue_key = (str(row[6]), int(row[7]), str(row[8]), segment)
        atom_names_by_residue.setdefault(residue_key, set()).add(str(row[3]))

    out_lines = []
    if ns.large:
        out_lines.append("LARGE FORMAT")
    out_lines.append("REMARK adjusted by smog_adjustPDB (python-native)")

    segment = 0
    for _, rec, serial, name, altloc, resn, chain, resi, icode, x, y, z, elem, ter_before, reset_before in recs:
        if ter_before:
            for _ in range(int(ter_before)):
                out_lines.append("TER")
                segment += 1
                cur_res = None
                if reset_before:
                    serial_out = 0
                    res_out = 0
        residue_key = (chain, resi, icode, segment)
        base = _base_residue_name(resn)
        if ns.map is None or ns.default:
            resn = _smog2_default_terminal_residue_name(
                resn,
                atom_names_by_residue.get(residue_key, set()),
                first_residue=residue_key in first_residues,
                last_residue=residue_key in last_residues,
            )
            base = _base_residue_name(resn)
        if residue_key in first_residues and base in rules.first_residue_names:
            resn = rules.first_residue_names[base]
        if residue_key in last_residues and base in rules.last_residue_names:
            resn = rules.last_residue_names[base]
        serial_out = serial if ns.PDBnums else serial_out + 1
        if ns.PDBnums:
            res_out = resi
        else:
            if cur_res != (chain, resi, icode, segment):
                res_out += 1
                cur_res = (chain, resi, icode, segment)
        out_lines.append(_fmt_atom(serial_out, name, altloc, resn, chain, res_out, " ", x, y, z, elem, rec))

    if ns.insertTER and pending_ter and out_lines and not out_lines[-1].startswith("TER"):
        out_lines.append("TER")

    Path(ns.o).write_text("\n".join(out_lines) + "\n")
    return 0
