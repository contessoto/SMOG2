from __future__ import annotations

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def _parse_ions_def(template: Path, ionnm: str):
    ionf = template / "ions.def"
    if not ionf.exists():
        matches = sorted(template.glob("*.ions.def"))
        if matches:
            ionf = matches[0]
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


def _section_name(line: str) -> str | None:
    s = line.strip()
    if s.startswith("[") and s.endswith("]"):
        return s.strip("[]").strip()
    return None


def _section_span(lines: list[str], name: str) -> tuple[int, int] | None:
    start = None
    for idx, line in enumerate(lines):
        section = _section_name(line)
        if section == name:
            start = idx
            continue
        if start is not None and section is not None:
            return start, idx
    return (start, len(lines)) if start is not None else None


def _fmt_top_number(value: float) -> str:
    return f"{float(value):g}"


def _fmt_xml_number(value: str) -> str:
    try:
        return f"{float(value):.5e}"
    except ValueError:
        return value


def _add_ion_to_top(top_in: Path, top_out: Path, ionnm: str, ionn: int, m: float, q: float, c12: float, c6: float) -> None:
    lines = top_in.read_text().splitlines(keepends=True)
    atomtypes = _section_span(lines, "atomtypes")
    if atomtypes is not None:
        _start, end = atomtypes
        has_ion = any(
            raw.strip()
            and not raw.lstrip().startswith(";")
            and raw.split()[0] == ionnm
            for raw in lines[atomtypes[0] + 1 : end]
        )
        if not has_ion:
            insert = end
            while insert > atomtypes[0] + 1 and not lines[insert - 1].strip():
                insert -= 1
            lines.insert(insert, f"{ionnm}\t{_fmt_top_number(m)}\t{_fmt_top_number(q)}\tA\t{_fmt_top_number(c6)}\t{_fmt_top_number(c12)}\n")

    first_moleculetype = next((idx for idx, raw in enumerate(lines) if _section_name(raw) == "moleculetype"), None)
    if first_moleculetype is not None:
        existing_text = "".join(lines)
        if f"\n{ionnm} 1\n" not in existing_text:
            block = [
                "[ moleculetype ]\n",
                f"{ionnm} 1\n",
                "\n",
                "[ atoms ]\n",
                f"1 {ionnm} 1 {ionnm} {ionnm} 1\n",
                "\n",
            ]
            lines[first_moleculetype:first_moleculetype] = block

    molecules = _section_span(lines, "molecules")
    if molecules is not None:
        _start, end = molecules
        if not any(raw.strip().split()[:1] == [ionnm] for raw in lines[molecules[0] + 1 : end] if raw.strip() and not raw.lstrip().startswith(";")):
            insert = end
            while insert > molecules[0] + 1 and not lines[insert - 1].strip():
                insert -= 1
            lines.insert(insert, f"{ionnm} {ionn}\n")

    top_out.write_text("".join(lines))


def _box_values(box: str) -> tuple[float, float, float]:
    values = [float(item) for item in box.split()[:3]]
    if len(values) != 3:
        raise SystemExit(f"Cannot parse GRO box line: {box}")
    return values[0], values[1], values[2]


def _ion_positions(box: str, count: int) -> list[tuple[float, float, float]]:
    dims = _box_values(box)
    ngrid = max(1, math.ceil(count ** (1 / 3)))
    spacing = tuple(dim / (ngrid + 1) for dim in dims)
    positions: list[tuple[float, float, float]] = []
    for ix in range(ngrid):
        for iy in range(ngrid):
            for iz in range(ngrid):
                positions.append((spacing[0] * (ix + 1), spacing[1] * (iy + 1), spacing[2] * (iz + 1)))
                if len(positions) == count:
                    return positions
    return positions


def _write_ion_gro(gro_in: Path, gro_out: Path, ionnm: str, ionn: int) -> None:
    head, atoms, box = _parse_gro(gro_in)
    start_idx = len(atoms) + 1
    for offset, (x, y, z) in enumerate(_ion_positions(box, ionn)):
        idx = start_idx + offset
        atoms.append(_fmt_gro_atom(idx, ionnm[:5], ionnm[:5], idx, x, y, z))
    out = [f"{head}. Note: ions added with smog_ions (version 2.4.5)", str(len(atoms)), *atoms, box]
    gro_out.write_text("\n".join(out) + "\n")


def _parse_extras_nonbond_rows(template: Path) -> list[tuple[str, str, list[str]]]:
    extras = next(iter(sorted(template.glob("*.extras"))), template / "extras")
    if not extras.exists():
        return []
    rows: list[tuple[str, str, list[str]]] = []
    for raw in extras.read_text().splitlines():
        text = raw.strip()
        if not text or text.startswith((";", "#")) or not text.startswith("nonbond_params"):
            continue
        cols = text.replace("<", " ").split()
        if len(cols) >= 6:
            rows.append((cols[1], cols[2], cols[4:]))
    return rows


def _pair_key(type1: str, type2: str) -> tuple[str, str]:
    return tuple(sorted((type1, type2)))


def _indent_xml(root: ET.Element) -> None:
    if hasattr(ET, "indent"):
        ET.indent(root, space=" ")


def _write_opensmog_xml_with_ion(xml_in: Path, xml_out: Path, template: Path | None, ionnm: str) -> None:
    root = ET.parse(xml_in).getroot()
    nonbond = root.find("nonbond")
    if nonbond is None:
        nonbond = ET.SubElement(root, "nonbond")
    bytype = nonbond.find("nonbond_bytype")
    if bytype is None:
        bytype = ET.SubElement(nonbond, "nonbond_bytype")
        ET.SubElement(bytype, "expression", {"expr": "null"})
        parameter = ET.SubElement(bytype, "parameter")
        parameter.text = "null"

    parameters = [node.text.strip() for node in bytype.findall("parameter") if node.text and node.text.strip()]
    if not parameters:
        parameters = ["null"]
        parameter = ET.SubElement(bytype, "parameter")
        parameter.text = "null"

    existing_params = [
        node for node in bytype.findall("nonbond_param")
        if node.attrib.get("type1") and node.attrib.get("type2")
    ]
    existing_types = {node.attrib["type1"] for node in existing_params} | {node.attrib["type2"] for node in existing_params}
    existing_types.add(ionnm)
    existing_pairs = {_pair_key(node.attrib["type1"], node.attrib["type2"]) for node in existing_params}

    template_rows = _parse_extras_nonbond_rows(template) if template is not None else []
    added_from_template = False
    for type1, type2, values in template_rows:
        if ionnm not in {type1, type2}:
            continue
        other = type2 if type1 == ionnm else type1
        if other not in existing_types:
            continue
        key = _pair_key(type1, type2)
        if key in existing_pairs:
            continue
        attrs = {"type1": type1, "type2": type2}
        attrs.update({parameter: _fmt_xml_number(value) for parameter, value in zip(parameters, values)})
        ET.SubElement(bytype, "nonbond_param", attrs)
        existing_pairs.add(key)
        existing_types.update((type1, type2))
        added_from_template = True

    if not added_from_template:
        for other in sorted(existing_types):
            key = _pair_key(other, ionnm)
            if key in existing_pairs:
                continue
            attrs = {"type1": other, "type2": ionnm}
            attrs.update({parameter: ("0" if parameter == "null" else "0.00000e+00") for parameter in parameters})
            ET.SubElement(bytype, "nonbond_param", attrs)
            existing_pairs.add(key)

    _indent_xml(root)
    xml_out.write_text(ET.tostring(root, encoding="unicode") + "\n")


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
    p.add_argument("-OpenSMOG", default=None)
    p.add_argument("-OpenSMOGout", default=None)
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_ions -f in.top -g in.gro -ionnm NAME -ionn N [-ionq Q -ionm M -ionC12 C12 -ionC6 C6 | -t template] [-OpenSMOG in.xml -OpenSMOGout out.xml]")
        return 1
    if ns.ionnm is None or ns.ionn < 1:
        raise SystemExit("Please indicate ion name with -ionnm and positive number with -ionn")

    if ns.t:
        m, q, c12, c6 = _parse_ions_def(Path(ns.t), ns.ionnm)
    else:
        if None in (ns.ionm, ns.ionq, ns.ionC12, ns.ionC6):
            raise SystemExit("Must provide -ionm -ionq -ionC12 -ionC6 when -t is not used")
        m, q, c12, c6 = ns.ionm, ns.ionq, ns.ionC12, ns.ionC6

    _add_ion_to_top(Path(ns.f), Path(ns.of), ns.ionnm, ns.ionn, m, q, c12, c6)
    _write_ion_gro(Path(ns.g), Path(ns.og), ns.ionnm, ns.ionn)
    if ns.OpenSMOG and ns.OpenSMOGout:
        _write_opensmog_xml_with_ion(Path(ns.OpenSMOG), Path(ns.OpenSMOGout), Path(ns.t) if ns.t else None, ns.ionnm)
    return 0


if __name__ == "__main__":
    import sys

    raise SystemExit(main(sys.argv[1:]))
