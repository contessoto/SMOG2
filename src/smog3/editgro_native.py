"""Native GRO box editing helper.

This module implements the tutorial-covered subset of ``smog_editgro``:
constructing buffered boxes, centering coordinates, and wrapping atoms through
periodic boundaries.  It operates directly on GRO text and does not call Perl or
SMOG2.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def _parse_gro(path: Path):
    lines = path.read_text().splitlines()
    head = lines[0]
    n = int(lines[1])
    atoms = lines[2 : 2 + n]
    box = [float(x) for x in lines[2 + n].split()]
    coords = []
    for ln in atoms:
        x = float(ln[20:28]); y = float(ln[28:36]); z = float(ln[36:44])
        coords.append((ln, x, y, z))
    return head, coords, box


def _write_gro(path: Path, head: str, coords, box):
    """Write GRO coordinates with updated positions and box dimensions."""

    out = [head, str(len(coords))]
    for ln, x, y, z in coords:
        out.append(f"{ln[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[44:]}")
    out.append(" ".join(f"{v:.5f}" for v in box))
    path.write_text("\n".join(out) + "\n")


def main(argv: list[str]) -> int:
    """Run native GRO centering, box construction, and PBC wrapping."""

    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-g", default="smog.gro")
    p.add_argument("-og", default="smog.box.gro")
    p.add_argument("-boxtype", default=None)
    p.add_argument("-d", type=float, default=None)
    p.add_argument("-c", action="store_true")
    p.add_argument("-pbc", action="store_true")
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_editgro -g in.gro -og out.gro [-boxtype cubic|octahedron|dodecahedron|rectangular] [-d dist] [-c] [-pbc]")
        return 1

    head, coords, box = _parse_gro(Path(ns.g))
    xs = [c[1] for c in coords]; ys = [c[2] for c in coords]; zs = [c[3] for c in coords]
    minx, maxx = min(xs), max(xs); miny, maxy = min(ys), max(ys); minz, maxz = min(zs), max(zs)
    if ns.boxtype:
        dist = 0.0 if ns.d is None else ns.d
        lx = (maxx - minx) + 2 * dist
        ly = (maxy - miny) + 2 * dist
        lz = (maxz - minz) + 2 * dist
        b = max(lx, ly, lz)
        bt = ns.boxtype.lower()
        if bt in {"cubic", "rectangular"}:
            box = [b, b, b]
        elif bt == "dodecahedron":
            box = [b, b, b, 0.0, 0.0, b / 2.0, 0.0, b / 2.0, 0.0]
        elif bt == "octahedron":
            box = [b, b, b, b / 3.0, b / 3.0, 0.0, b / 3.0, 0.0, b / 3.0]
        else:
            raise SystemExit(f"boxtype value \"{ns.boxtype}\" is not supported.")

    if ns.c:
        cx = (minx + maxx) / 2.0; cy = (miny + maxy) / 2.0; cz = (minz + maxz) / 2.0
        tx, ty, tz = box[0] / 2.0, box[1] / 2.0, box[2] / 2.0
        dx, dy, dz = tx - cx, ty - cy, tz - cz
        coords = [(ln, x + dx, y + dy, z + dz) for (ln, x, y, z) in coords]

    if ns.pbc:
        lx, ly, lz = box[0], box[1], box[2]
        wrapped = []
        for ln, x, y, z in coords:
            while x < 0: x += lx
            while y < 0: y += ly
            while z < 0: z += lz
            while x >= lx: x -= lx
            while y >= ly: y -= ly
            while z >= lz: z -= lz
            wrapped.append((ln, x, y, z))
        coords = wrapped

    _write_gro(Path(ns.og), head, coords, box)
    return 0
