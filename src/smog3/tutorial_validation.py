"""Public SMOG tutorial validation runner.

This validation-only module reads ``validation/tutorials/tutorial_manifest.yml``
and executes tutorial model-generation workflows on two sides: official SMOG2
inside Docker for the baseline, and local or installed SMOG3 for the candidate.
It compares generated topology, coordinate, index, contact, and OpenSMOG XML
files using the same documented comparator policy as SMOG-CHECK parity.
"""

from __future__ import annotations

import argparse
import difflib
import itertools
import json
import os
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .parity_direct import _compare_file, _drop_top_header_metadata, compare_existing_dirs

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "validation" / "tutorials" / "tutorial_manifest.yml"
DEFAULT_RUN_ROOT = ROOT / "validation" / "tutorials" / "runs"
DOCKER_IMAGE = "smogserver/smog2:stable"
COMPARATOR_POLICY = (
    "Exact output comparison with the documented SMOG-CHECK parity policy: generated "
    "topology/XML metadata and tiny topology floating-point print noise may be ignored, "
    "but scientific section values are not relaxed."
)
BAD_STATUSES = {"DIFF", "SMOG2_ERROR", "SMOG3_ERROR", "UNSUPPORTED_BY_SMOG3"}
ALL_STATUSES = (
    "PASS",
    "DIFF",
    "SMOG2_ERROR",
    "SMOG3_ERROR",
    "MISSING_INPUT",
    "MISSING_DOWNLOAD",
    "MANUAL_INPUT_REQUIRED",
    "UNSUPPORTED_BY_SMOG3",
    "NOT_GENERATION_TEST",
)


@dataclass(frozen=True)
class TutorialCase:
    """A single public-tutorial model-generation validation case."""

    raw: dict[str, Any]

    @property
    def case_id(self) -> str:
        return str(self.raw["id"])

    @property
    def status(self) -> str:
        return str(self.raw.get("status", "missing_input"))

    @property
    def implemented(self) -> bool:
        return self.status == "implemented"

    @property
    def workflow(self) -> bool:
        return bool(self.raw.get("workflow_commands"))

    @property
    def expected_outputs(self) -> list[str]:
        return [str(item) for item in self.raw.get("expected_outputs", [])]

    @property
    def include_xml(self) -> bool:
        return "model.xml" in self.expected_outputs

    @property
    def coord_name(self) -> str:
        return "model.g96" if "model.g96" in self.expected_outputs else "model.gro"


def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    """Load the tutorial manifest.

    The file is intentionally a YAML-compatible JSON subset. Using JSON keeps the
    runner dependency-free while still giving YAML tooling a valid manifest.
    """

    return json.loads(path.read_text(encoding="utf-8"))


def tutorial_cases(path: Path = DEFAULT_MANIFEST) -> list[TutorialCase]:
    """Return all tutorial manifest cases in manifest order."""

    return [TutorialCase(raw) for raw in load_manifest(path).get("cases", [])]


def _rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _shell_join(args: list[str]) -> str:
    return " ".join(subprocess.list2cmdline([arg]) for arg in args)


def _read_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def _first_diff(a: Path, b: Path, limit: int = 80) -> str:
    if not a.exists() or not b.exists():
        return ""
    left_lines = _read_lines(a)
    right_lines = _read_lines(b)
    if max(len(left_lines), len(right_lines)) > 50000:
        for idx, (left, right) in enumerate(itertools.zip_longest(left_lines, right_lines, fillvalue="<missing>"), start=1):
            if left != right:
                return f"large-file first differing line {idx}\n--- {a}\n+++ {b}\n- {left}\n+ {right}"
        return "large files differ"
    diff = difflib.unified_diff(left_lines, right_lines, fromfile=str(a), tofile=str(b), n=3, lineterm="")
    return "\n".join(itertools.islice(diff, limit))


def _topology_section_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    section: str | None = None
    for line in _read_lines(path):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped.strip("[]").strip()
            counts.setdefault(section, 0)
        elif section and stripped and not stripped.startswith(";"):
            counts[section] += 1
    return counts


def _count_contacts(path: Path) -> int | None:
    if not path.exists():
        return None
    return sum(1 for line in _read_lines(path) if line.strip() and not line.lstrip().startswith(";"))


def _coord_atom_count(path: Path) -> int | None:
    lines = _read_lines(path)
    if path.suffix == ".gro" and len(lines) >= 2:
        try:
            return int(lines[1].strip())
        except ValueError:
            return None
    if path.suffix == ".g96":
        in_position = False
        count = 0
        for line in lines:
            stripped = line.strip()
            if stripped == "POSITION":
                in_position = True
                continue
            if in_position and stripped == "END":
                return count
            if in_position and stripped:
                count += 1
        return count or None
    return None


def _xml_status(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False, "parse_ok": False}
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        return {"exists": True, "parse_ok": False, "error": str(exc)}
    return {"exists": True, "parse_ok": True, "root": root.tag}


def _gro_atom_records(path: Path) -> tuple[str, list[str], str] | None:
    lines = _read_lines(path)
    if len(lines) < 3:
        return None
    try:
        count = int(lines[1].strip())
    except ValueError:
        return None
    if len(lines) < count + 3:
        return None
    return lines[0], lines[2 : 2 + count], lines[2 + count]


def _gro_residue_name(line: str) -> str:
    return line[5:10].strip()


def _gro_coordinates(line: str) -> tuple[float, float, float] | None:
    try:
        return float(line[20:28]), float(line[28:36]), float(line[36:44])
    except ValueError:
        return None


def _compare_ion_gro(left: Path, right: Path) -> dict[str, Any]:
    left_records = _gro_atom_records(left)
    right_records = _gro_atom_records(right)
    if left_records is None or right_records is None:
        return _compare_file(left, right)
    _left_header, left_atoms, left_box = left_records
    _right_header, right_atoms, right_box = right_records
    if len(left_atoms) != len(right_atoms):
        return {"match": False, "reason": "GRO atom count differs"}
    if left_box != right_box:
        return {"match": False, "reason": "GRO box line differs"}

    ion_names = {"CL", "K", "MG", "NA"}
    left_non_ions = [line for line in left_atoms if _gro_residue_name(line) not in ion_names]
    right_non_ions = [line for line in right_atoms if _gro_residue_name(line) not in ion_names]
    if left_non_ions != right_non_ions:
        return {"match": False, "reason": "non-ion GRO atom records differ", "diff": _first_diff(left, right)}

    left_ions = Counter(_gro_residue_name(line) for line in left_atoms if _gro_residue_name(line) in ion_names)
    right_ions = Counter(_gro_residue_name(line) for line in right_atoms if _gro_residue_name(line) in ion_names)
    if left_ions != right_ions:
        return {"match": False, "reason": "ion species/counts differ", "baseline_ions": dict(left_ions), "candidate_ions": dict(right_ions)}

    box = [float(value) for value in right_box.split()[:3]]
    for line in right_atoms:
        if _gro_residue_name(line) not in ion_names:
            continue
        coords = _gro_coordinates(line)
        if coords is None or any(not (0.0 <= coord <= box[idx]) for idx, coord in enumerate(coords)):
            return {"match": False, "reason": "candidate ion coordinate outside finite box", "line": line}
    return {
        "match": True,
        "ignored": "smog_ions places ions stochastically; compared non-ion GRO records exactly, ion species/counts exactly, and candidate ion coordinates inside the same box",
        "ion_counts": dict(right_ions),
    }


def _xml_signature(node: ET.Element) -> tuple[Any, ...]:
    text = (node.text or "").strip()
    children = [_xml_signature(child) for child in list(node)]
    if node.tag == "nonbond_bytype":
        fixed = [child for child in children if child[0] != "nonbond_param"]
        params = sorted(child for child in children if child[0] == "nonbond_param")
        children = fixed + params
    return (node.tag, tuple(sorted(node.attrib.items())), text, tuple(children))


def _compare_xml_semantic(left: Path, right: Path) -> dict[str, Any]:
    if not left.exists() or not right.exists():
        return _compare_file(left, right)
    try:
        left_sig = _xml_signature(ET.parse(left).getroot())
        right_sig = _xml_signature(ET.parse(right).getroot())
    except ET.ParseError as exc:
        return {"match": False, "reason": f"XML parse error: {exc}"}
    if left_sig == right_sig:
        return {"match": True, "ignored": "OpenSMOG XML comments/whitespace and nonbond_param order"}
    return {"match": False, "reason": "OpenSMOG XML semantic content differs", "diff": _first_diff(left, right)}


def _normalize_atomtype_layout(lines: list[str]) -> list[str]:
    normalized: list[str] = []
    section: str | None = None
    for line in lines:
        name = line.strip()
        if name.startswith("[") and name.endswith("]"):
            section = name.strip("[]").strip()
            normalized.append(f"[ {section} ]")
            continue
        if section == "atomtypes" and line.strip() and not line.lstrip().startswith(";"):
            normalized.append("\t".join(line.split()))
        else:
            normalized.append(line)
    return normalized


def _compare_ion_top(left: Path, right: Path) -> dict[str, Any]:
    standard = _compare_file(left, right)
    if standard.get("match"):
        return standard
    if not left.exists() or not right.exists():
        return standard
    left_norm = _normalize_atomtype_layout(_read_lines(left))
    right_norm = _normalize_atomtype_layout(_read_lines(right))
    if _normalize_atomtype_layout(_drop_top_header_metadata(left_norm)) == _normalize_atomtype_layout(_drop_top_header_metadata(right_norm)):
        return {"match": True, "ignored": "topology header metadata and atomtypes whitespace rewritten by smog_ions"}
    return standard


def _file_diagnostics(baseline_dir: Path, candidate_dir: Path, outputs: list[str]) -> dict[str, Any]:
    diagnostics: dict[str, Any] = {}
    for name in outputs:
        left = baseline_dir / name
        right = candidate_dir / name
        item: dict[str, Any] = {
            "baseline_exists": left.exists(),
            "candidate_exists": right.exists(),
            "baseline_lines": len(_read_lines(left)) if left.exists() else None,
            "candidate_lines": len(_read_lines(right)) if right.exists() else None,
        }
        if left.exists() and right.exists() and left.read_bytes() != right.read_bytes():
            item["first_diff"] = _first_diff(left, right)
        if Path(name).suffix in {".gro", ".g96"}:
            item["baseline_atom_count"] = _coord_atom_count(left)
            item["candidate_atom_count"] = _coord_atom_count(right)
        if Path(name).suffix == ".contacts":
            item["baseline_contact_lines"] = _count_contacts(left)
            item["candidate_contact_lines"] = _count_contacts(right)
        if Path(name).suffix == ".top":
            item["baseline_section_counts"] = _topology_section_counts(left)
            item["candidate_section_counts"] = _topology_section_counts(right)
        if Path(name).suffix == ".xml":
            item["baseline_xml"] = _xml_status(left)
            item["candidate_xml"] = _xml_status(right)
        diagnostics[name] = item
    return diagnostics


def _compare_outputs(baseline_dir: Path, candidate_dir: Path, outputs: list[str], case: TutorialCase | None = None) -> dict[str, Any]:
    standard = {"model.top", "model.gro", "model.ndx", "model.contacts", "model.xml"}
    ion_workflow = bool(case and case.raw.get("ion_workflow"))
    if set(outputs).issubset(standard) and not ion_workflow:
        report = compare_existing_dirs(baseline_dir, candidate_dir, include_xml="model.xml" in outputs)
        report["comparisons"] = {name: report["comparisons"][name] for name in outputs}
        report["ok"] = all(item.get("match") for item in report["comparisons"].values())
        return report

    comparisons = {}
    for name in outputs:
        left = baseline_dir / name
        right = candidate_dir / name
        if ion_workflow and Path(name).suffix == ".gro":
            comparisons[name] = _compare_ion_gro(left, right)
        elif ion_workflow and Path(name).suffix == ".xml":
            comparisons[name] = _compare_xml_semantic(left, right)
        elif ion_workflow and Path(name).suffix == ".top":
            comparisons[name] = _compare_ion_top(left, right)
        else:
            comparisons[name] = _compare_file(left, right)
    return {
        "baseline_dir": str(baseline_dir),
        "candidate_dir": str(candidate_dir),
        "comparisons": comparisons,
        "ok": all(item.get("match") for item in comparisons.values()),
    }


def _copy_inputs(case: TutorialCase, run_inputs: Path) -> tuple[dict[str, str], list[str]]:
    missing: list[str] = []
    context: dict[str, str] = {}
    case_input_dir = run_inputs / case.case_id
    case_input_dir.mkdir(parents=True, exist_ok=True)

    for rel in case.raw.get("required_input_files", []):
        src = ROOT / str(rel)
        if not src.exists():
            missing.append(str(rel))
            continue
        dst = case_input_dir / src.name
        shutil.copy2(src, dst)
        if str(rel) == str(case.raw.get("input_pdb")):
            context["input_pdb_rel"] = _rel(dst)
            context["input_pdb_docker"] = f"/workdir/{_rel(dst)}"
        if str(rel) == str(case.raw.get("user_contacts")):
            context["user_contacts_rel"] = _rel(dst)
            context["user_contacts_docker"] = f"/workdir/{_rel(dst)}"

    return context, missing


def _copy_workflow_inputs(case: TutorialCase, work_dir: Path) -> list[str]:
    """Copy downloaded tutorial assets into an isolated workflow directory."""

    missing: list[str] = []
    mappings = case.raw.get("workflow_input_files") or [
        {"source": rel, "dest": Path(str(rel)).name}
        for rel in case.raw.get("required_input_files", [])
    ]
    for item in mappings:
        source = ROOT / str(item["source"])
        if not source.exists():
            missing.append(str(item["source"]))
            continue
        dest = work_dir / str(item["dest"])
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
    return missing


def _expand_args(raw_args: list[str], context: dict[str, str], docker: bool) -> list[str]:
    values = {
        "user_contacts": context.get("user_contacts_docker" if docker else "user_contacts_rel", ""),
        "model_xml": context["model_xml_docker" if docker else "model_xml_rel"],
    }
    return [str(arg).format(**values) for arg in raw_args]


def _base_smog2_args(case: TutorialCase, context: dict[str, str], baseline_dir: Path) -> list[str]:
    context = dict(context)
    context["model_xml_docker"] = f"/workdir/{_rel(baseline_dir / 'model.xml')}"
    return [
        "smog2",
        "-i",
        context["input_pdb_docker"],
        "-SCMorig",
        "-keep4SCM",
        "-o",
        f"/workdir/{_rel(baseline_dir / 'model.top')}",
        "-g",
        f"/workdir/{_rel(baseline_dir / case.coord_name)}",
        "-n",
        f"/workdir/{_rel(baseline_dir / 'model.ndx')}",
        "-s",
        f"/workdir/{_rel(baseline_dir / 'model.contacts')}",
        *_expand_args([str(arg) for arg in case.raw.get("smog2_args", [])], context, docker=True),
    ]


def _base_smog3_args(case: TutorialCase, context: dict[str, str], candidate_dir: Path, use_installed: bool) -> list[str]:
    context = dict(context)
    context["model_xml_rel"] = _rel(candidate_dir / "model.xml")
    command = ["smog3"] if use_installed else ["python3", "-m", "smog3.smog2_native"]
    return [
        *command,
        "-i",
        context["input_pdb_rel"],
        "-o",
        _rel(candidate_dir / "model.top"),
        "-g",
        _rel(candidate_dir / case.coord_name),
        "-n",
        _rel(candidate_dir / "model.ndx"),
        "-s",
        _rel(candidate_dir / "model.contacts"),
        *_expand_args([str(arg) for arg in case.raw.get("smog3_args", [])], context, docker=False),
    ]


def _run_smog2(case: TutorialCase, context: dict[str, str], baseline_dir: Path, image: str) -> tuple[int, str, list[str]]:
    args = _base_smog2_args(case, context, baseline_dir)
    script = f"set -euo pipefail\ncd /workdir\n{_shell_join(args)}\n"
    proc = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{ROOT}:/workdir", image, "bash", "-lc", script],
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr, args


def _run_smog3(
    case: TutorialCase,
    context: dict[str, str],
    candidate_dir: Path,
    use_installed: bool,
    no_perl_bin: Path,
    perl_log: Path,
) -> tuple[int, str, list[str]]:
    args = _base_smog3_args(case, context, candidate_dir, use_installed)
    env = os.environ.copy()
    env["PATH"] = f"{no_perl_bin}:{env.get('PATH', '')}"
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    env["SMOG3_USE_SCM_DEFAULTS"] = "1"
    env["SMOG3_PERL_SENTINEL_LOG"] = str(perl_log)
    if not use_installed:
        env["PYTHONPATH"] = f"{ROOT / 'src'}:{env.get('PYTHONPATH', '')}"
    proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout + proc.stderr, args


def _workflow_adjust_command() -> str:
    return (
        "python3 -c 'from smog3.adjustpdb_native import main; "
        "import sys; raise SystemExit(main(sys.argv[1:]))'"
    )


def _workflow_smog3_command(_use_installed: bool) -> str:
    return "python3 -m smog3.smogcheck_dropin_smog2"


def _workflow_ions_command(use_installed: bool) -> str:
    return "smog3-ions" if use_installed else "python3 -m smog3.ions_native"


def _workflow_extract_command(use_installed: bool) -> str:
    return "smog3-extract" if use_installed else "python3 -m smog3.extract_native"


def _workflow_candidate_command(command: str, use_installed: bool) -> str:
    if "smog_adjustPDB " in command:
        return command.replace("smog_adjustPDB", _workflow_adjust_command(), 1)
    if "smog_extract " in command:
        return command.replace("smog_extract", _workflow_extract_command(use_installed), 1)
    if "smog_ions " in command:
        return command.replace("smog_ions", _workflow_ions_command(use_installed), 1)
    if "smog2 " in command:
        return command.replace("smog2", _workflow_smog3_command(use_installed), 1)
    return command


def _run_workflow_smog2(case: TutorialCase, baseline_dir: Path, image: str) -> tuple[int, str, list[str]]:
    commands = [str(cmd) for cmd in case.raw.get("workflow_commands", [])]
    script = "set -euo pipefail\n" f"cd /workdir/{_rel(baseline_dir)}\n" + "\n".join(commands) + "\n"
    proc = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{ROOT}:/workdir", image, "bash", "-lc", script],
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr, commands


def _run_workflow_smog3(
    case: TutorialCase,
    candidate_dir: Path,
    use_installed: bool,
    no_perl_bin: Path,
    perl_log: Path,
) -> tuple[int, str, list[str]]:
    commands = [_workflow_candidate_command(str(cmd), use_installed) for cmd in case.raw.get("workflow_commands", [])]
    env = os.environ.copy()
    env["PATH"] = f"{no_perl_bin}:{env.get('PATH', '')}"
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    env["SMOG3_USE_SCM_DEFAULTS"] = "1"
    env["SMOG3_PERL_SENTINEL_LOG"] = str(perl_log)
    if not use_installed:
        env["PYTHONPATH"] = f"{ROOT / 'src'}:{env.get('PYTHONPATH', '')}"
    proc = subprocess.run(
        ["bash", "-lc", "set -euo pipefail\n" + "\n".join(commands) + "\n"],
        cwd=candidate_dir,
        capture_output=True,
        text=True,
        env=env,
    )
    return proc.returncode, proc.stdout + proc.stderr, commands


def _write_no_perl_sentinel(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    sentinel = path / "perl"
    sentinel.write_text(
        "#!/usr/bin/env bash\n"
        "echo \"SMOG3 attempted to invoke perl: $*\" >> \"${SMOG3_PERL_SENTINEL_LOG:-/tmp/smog3-tutorial-perl.log}\"\n"
        "exit 127\n",
        encoding="utf-8",
    )
    sentinel.chmod(0o755)


def _preclassified_report(case: TutorialCase, status: str, reports_dir: Path, reason: str = "") -> dict[str, Any]:
    report = {
        "case": case.case_id,
        "tutorial_name": case.raw.get("tutorial_name"),
        "status": status,
        "model_type": case.raw.get("model_type"),
        "source_url": case.raw.get("source_url"),
        "required_input_files": case.raw.get("required_input_files", []),
        "reason": reason or case.raw.get("notes", ""),
        "comparison_policy": case.raw.get("comparison_policy", COMPARATOR_POLICY),
        "smog2_command": None,
        "smog3_command": None,
        "comparison": None,
        "diagnostics": {},
    }
    (reports_dir / f"{case.case_id}.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def _missing_status(missing: list[str]) -> str:
    for item in missing:
        if item.startswith("validation/tutorials/assets/"):
            return "MISSING_DOWNLOAD"
    return "MISSING_INPUT"


def _generation_outputs_exist(directory: Path, outputs: list[str]) -> bool:
    # SMOG2 consumes user-provided contact files into the topology and does not
    # always emit model.contacts.  For one-sided generation checks, require the
    # structural outputs and treat model.contacts as optional.
    required = [name for name in outputs if name != "model.contacts"]
    return all((directory / name).exists() for name in required)


def _run_case(
    case: TutorialCase,
    run_root: Path,
    image: str,
    use_installed: bool,
    with_smog2_baseline: bool,
    smog3_only: bool,
    smog2_only: bool,
    no_perl_bin: Path,
    perl_log: Path,
) -> dict[str, Any]:
    reports_dir = run_root / "reports"
    if not case.implemented:
        return _preclassified_report(case, case.status.upper(), reports_dir)

    baseline_dir = run_root / "smog2_baseline" / case.case_id
    candidate_dir = run_root / "smog3_candidate" / case.case_id
    baseline_dir.mkdir(parents=True, exist_ok=True)
    candidate_dir.mkdir(parents=True, exist_ok=True)
    if case.workflow:
        missing = _copy_workflow_inputs(case, baseline_dir) + _copy_workflow_inputs(case, candidate_dir)
        context: dict[str, str] = {}
    else:
        context, missing = _copy_inputs(case, run_root / "inputs")
    if missing:
        return _preclassified_report(case, _missing_status(missing), reports_dir, reason="Missing files: " + ", ".join(missing))

    brc = 0
    bout = ""
    bargs: list[str] | None = None
    logs_dir = run_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    if with_smog2_baseline and not smog3_only:
        if case.workflow:
            brc, bout, bargs = _run_workflow_smog2(case, baseline_dir, image)
        else:
            brc, bout, bargs = _run_smog2(case, context, baseline_dir, image)
        (baseline_dir / "smog2.log").write_text(bout, encoding="utf-8", errors="replace")
        (logs_dir / f"{case.case_id}.smog2.log").write_text(bout, encoding="utf-8", errors="replace")
    else:
        bargs = None

    crc: int | None = None
    cout = ""
    cargs: list[str] | None = None
    if not smog2_only:
        if case.workflow:
            crc, cout, cargs = _run_workflow_smog3(case, candidate_dir, use_installed, no_perl_bin, perl_log)
        else:
            crc, cout, cargs = _run_smog3(case, context, candidate_dir, use_installed, no_perl_bin, perl_log)
        (candidate_dir / "smog3.log").write_text(cout, encoding="utf-8", errors="replace")
        (logs_dir / f"{case.case_id}.smog3.log").write_text(cout, encoding="utf-8", errors="replace")

    if brc != 0:
        status = "SMOG2_ERROR"
        comparison = _compare_outputs(baseline_dir, candidate_dir, case.expected_outputs, case)
    elif smog2_only:
        generated = _generation_outputs_exist(baseline_dir, case.expected_outputs)
        status = "PASS" if generated else "SMOG2_ERROR"
        comparison = {"ok": generated, "mode": "smog2_only_generation_check", "comparisons": {}}
    elif crc != 0:
        status = "SMOG3_ERROR"
        comparison = _compare_outputs(baseline_dir, candidate_dir, case.expected_outputs, case)
    elif smog3_only or not with_smog2_baseline:
        generated = _generation_outputs_exist(candidate_dir, case.expected_outputs)
        status = "PASS" if generated else "SMOG3_ERROR"
        comparison = {"ok": generated, "mode": "smog3_only_generation_check", "comparisons": {}}
    else:
        comparison = _compare_outputs(baseline_dir, candidate_dir, case.expected_outputs, case)
        status = "PASS" if comparison.get("ok") else "DIFF"

    report = {
        "case": case.case_id,
        "tutorial_name": case.raw.get("tutorial_name"),
        "status": status,
        "model_type": case.raw.get("model_type"),
        "source_url": case.raw.get("source_url"),
        "input_source": case.raw.get("input_source"),
        "required_input_files": case.raw.get("required_input_files", []),
        "expected_outputs": case.expected_outputs,
        "baseline_dir": str(baseline_dir),
        "candidate_dir": str(candidate_dir),
        "baseline_returncode": brc,
        "candidate_returncode": crc,
        "baseline_tail": bout[-2000:],
        "candidate_tail": cout[-2000:],
        "smog2_command": bargs,
        "smog3_command": cargs,
        "comparison_policy": case.raw.get("comparison_policy", COMPARATOR_POLICY),
        "comparison": comparison,
        "diagnostics": _file_diagnostics(baseline_dir, candidate_dir, case.expected_outputs),
    }
    (reports_dir / f"{case.case_id}.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def _smog2_version(image: str) -> str:
    proc = subprocess.run(
        ["docker", "run", "--rm", image, "bash", "-lc", "smog2 -v 2>&1 | head -20"],
        capture_output=True,
        text=True,
    )
    return " ".join((proc.stdout + proc.stderr).split()) or "unknown"


def _render_summary(summary: dict[str, Any]) -> str:
    lines = [
        "# SMOG3 Tutorial Validation Summary",
        "",
        f"- Run directory: `{summary['run_root']}`",
        f"- SMOG3 version: `{summary['smog3_version']}`",
        f"- SMOG2 Docker image: `{summary['docker_image']}`",
        f"- SMOG2 Docker version: `{summary['smog2_version']}`",
        f"- SMOG3 Perl invocations: `{summary['smog3_perl_invocations']}`",
        f"- Comparator policy: {COMPARATOR_POLICY}",
        "",
        "| Status | Count |",
        "| --- | ---: |",
    ]
    for status, count in summary["counts"].items():
        lines.append(f"| {status} | {count} |")
    lines.extend(["", "| Case | Status | Tutorial Category | Report |", "| --- | --- | --- | --- |"])
    for case in summary["cases"]:
        lines.append(f"| {case['case']} | {case['status']} | {case.get('tutorial_name', '')} | `{case['report']}` |")
    lines.append("")
    return "\n".join(lines)


def _write_summary(run_root: Path, reports: list[dict[str, Any]], image: str, perl_log: Path) -> dict[str, Any]:
    counter = Counter(report.get("status", "UNKNOWN") for report in reports)
    counts = {status: counter.get(status, 0) for status in ALL_STATUSES}
    perl_count = sum(1 for _line in perl_log.open(encoding="utf-8")) if perl_log.exists() else 0
    summary = {
        "run_root": str(run_root),
        "case_count": len(reports),
        "counts": counts,
        "smog3_version": __version__,
        "docker_image": image,
        "smog2_version": _smog2_version(image),
        "smog3_perl_invocations": perl_count,
        "smog3_perl_log": str(perl_log),
        "comparator_policy": COMPARATOR_POLICY,
        "cases": [
            {
                "case": report["case"],
                "status": report["status"],
                "tutorial_name": report.get("tutorial_name"),
                "report": str(run_root / "reports" / f"{report['case']}.json"),
            }
            for report in reports
        ],
    }
    summary_json = json.dumps(summary, indent=2) + "\n"
    summary_md = _render_summary(summary)
    (run_root / "tutorial_validation_summary.json").write_text(summary_json, encoding="utf-8")
    (run_root / "tutorial_validation_summary.md").write_text(summary_md, encoding="utf-8")
    (run_root / "tutorial_compare_summary.json").write_text(summary_json, encoding="utf-8")
    (run_root / "tutorial_compare_summary.md").write_text(summary_md, encoding="utf-8")
    return summary


def _select_cases(all_cases: list[TutorialCase], ns: argparse.Namespace) -> list[TutorialCase]:
    if ns.case:
        wanted = {item.strip() for value in ns.case for item in value.split(",") if item.strip()}
        selected = [case for case in all_cases if case.case_id in wanted]
        missing = sorted(wanted - {case.case_id for case in selected})
        if missing:
            raise SystemExit(f"Unknown tutorial case id(s): {', '.join(missing)}")
        return selected
    if ns.all:
        return all_cases
    return [case for case in all_cases if case.implemented]


def _new_run_root(base: Path) -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S")
    run_root = base / stamp
    suffix = 0
    while run_root.exists():
        suffix += 1
        run_root = base / f"{stamp}.{suffix}"
    return run_root


def run_validation(ns: argparse.Namespace) -> int:
    """Run all requested tutorial cases and write JSON/Markdown summaries."""

    if ns.download_first:
        fetch_script = ROOT / "scripts" / "fetch_smog_tutorial_assets.py"
        subprocess.run([sys.executable, str(fetch_script)], cwd=ROOT, check=True)
    cases = _select_cases(tutorial_cases(Path(ns.manifest)), ns)
    run_root = Path(ns.out_dir) if ns.out_dir else _new_run_root(Path(ns.run_root))
    for subdir in ("inputs", "smog2_baseline", "smog3_candidate", "logs", "reports"):
        (run_root / subdir).mkdir(parents=True, exist_ok=True)
    no_perl_bin = run_root / "no-perl-bin"
    perl_log = run_root / "smog3-perl-invocations.log"
    _write_no_perl_sentinel(no_perl_bin)
    perl_log.write_text("", encoding="utf-8")

    reports: list[dict[str, Any]] = []
    for case in cases:
        print(f"[tutorial] {case.case_id}: {case.raw.get('tutorial_name')} ({case.status})", flush=True)
        report = _run_case(
            case,
            run_root,
            ns.image,
            ns.use_installed_smog3,
            ns.with_smog2_baseline,
            ns.smog3_only,
            ns.smog2_only,
            no_perl_bin,
            perl_log,
        )
        print(f"[tutorial] {case.case_id}: {report['status']}", flush=True)
        reports.append(report)
        if not ns.keep_going and report["status"] in BAD_STATUSES:
            break

    summary = _write_summary(run_root, reports, ns.image, perl_log)
    print(_render_summary(summary))
    if summary["smog3_perl_invocations"]:
        return 3
    return 1 if any(report["status"] in BAD_STATUSES for report in reports) else 0


def list_cases(ns: argparse.Namespace) -> int:
    """Print the tutorial manifest cases with their automation status."""

    for case in tutorial_cases(Path(ns.manifest)):
        print(f"{case.case_id}\t{case.status}\t{case.raw.get('tutorial_name')}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for tutorial listing and comparison runs."""

    parser = argparse.ArgumentParser(description="Validate SMOG3 against public SMOG tutorial model-generation examples.")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--run-root", default=str(DEFAULT_RUN_ROOT))
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--image", default=DOCKER_IMAGE)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--use-installed-smog3", action="store_true")
    parser.add_argument("--with-smog2-baseline", action="store_true", default=True)
    parser.add_argument("--smog3-only", action="store_true")
    parser.add_argument("--smog2-only", action="store_true")
    parser.add_argument("--download-first", action="store_true")
    parser.add_argument("--keep-going", action="store_true", default=True)
    ns = parser.parse_args(argv)
    if ns.list:
        return list_cases(ns)
    return run_validation(ns)


if __name__ == "__main__":
    raise SystemExit(main())
