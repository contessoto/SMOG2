"""Small GROMACS-format parsing helpers used by SMOG3 tools.

These utilities preserve section text rather than building a full molecular
model.  They are shared by runtime helper tools such as extraction and energy
scaling, where SMOG2-compatible ordering and comments matter for parity.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class TopSection:
    """A named GROMACS topology section plus its raw body lines."""

    name: str
    lines: list[str]


def parse_ndx(path: str | Path) -> dict[str, list[int]]:
    """Parse an NDX file into group-name to atom-index lists."""

    groups: dict[str, list[int]] = {}
    current: str | None = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            groups[current] = []
            continue
        if current is None:
            continue
        groups[current].extend(int(tok) for tok in line.split())
    return groups


def parse_top_sections(path: str | Path) -> list[TopSection]:
    """Split a topology into preamble and bracketed sections without reformatting."""

    sections: list[TopSection] = []
    current = TopSection(name="__preamble__", lines=[])
    sections.append(current)
    for raw in Path(path).read_text().splitlines(keepends=True):
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current = TopSection(name=stripped[1:-1].strip(), lines=[])
            sections.append(current)
        else:
            current.lines.append(raw)
    return sections


def write_top_sections(path: str | Path, sections: list[TopSection]) -> None:
    """Write topology sections back in order, preserving raw section bodies."""

    out: list[str] = []
    for i, sec in enumerate(sections):
        if i > 0:
            out.append(f"[ {sec.name} ]\n")
        out.extend(sec.lines)
    Path(path).write_text("".join(out))
