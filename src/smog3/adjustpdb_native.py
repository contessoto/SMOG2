from __future__ import annotations

import argparse
from pathlib import Path


def _parse_map(path: Path) -> dict[tuple[str, str], tuple[str, str]]:
    mapping = {}
    if not path.exists():
        return mapping
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith(";") or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) >= 4:
            mapping[(parts[0], parts[1])] = (parts[2], parts[3])
    return mapping


def _fmt_atom(serial: int, name: str, altloc: str, resn: str, chain: str, resi: int, x: float, y: float, z: float, elem: str, rec: str="ATOM") -> str:
    return f"{rec:<6}{serial:>5d} {name:<4}{altloc:1}{resn:>3} {chain:1}{resi:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {elem:>2}"


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
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra or not ns.i:
        print("usage: smog_adjustPDB -i input.pdb [-o adjusted.pdb] [-map file] [-altloc] [-removeH] [-removewater] [-PDBnums] [-sort] [-large]")
        return 1

    m = _parse_map(Path(ns.map)) if ns.map else {}

    recs = []
    serial_out = 0
    cur_res = None
    res_out = 0
    lines = Path(ns.i).read_text().splitlines()
    for ln in lines:
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

        key = (resn, name)
        if key in m:
            resn, name = m[key]

        if ns.sort:
            sort_key = (chain, resi, name, serial)
        else:
            sort_key = (len(recs),)

        recs.append((sort_key, rec, serial, name, altloc if ns.altloc else " ", resn, chain, resi, x, y, z, elem))

    recs.sort(key=lambda t: t[0])
    out_lines = []
    if ns.large:
        out_lines.append("LARGE FORMAT")
    out_lines.append("REMARK adjusted by smog_adjustPDB (python-native)")

    for _, rec, serial, name, altloc, resn, chain, resi, x, y, z, elem in recs:
        serial_out = serial if ns.PDBnums else serial_out + 1
        if ns.PDBnums:
            res_out = resi
        else:
            if cur_res != (chain, resi):
                res_out += 1
                cur_res = (chain, resi)
        out_lines.append(_fmt_atom(serial_out, name, altloc, resn, chain, res_out, x, y, z, elem, rec))

    Path(ns.o).write_text("\n".join(out_lines) + "\n")
    return 0
