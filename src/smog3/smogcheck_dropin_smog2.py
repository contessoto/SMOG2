from __future__ import annotations

import json
import os
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from .smog2_native import main as smog2_native_main


MODEL_FLAGS = {
    "-AA",
    "-CA",
    "-AA2cg",
    "-AAnbcr2",
    "-AAgaussian",
    "-CAgaussian",
    "-AACC1",
    "-AACCD",
    "-AADIHE",
    "-AADIHE4",
    "-AAMATCH",
    "-AABOND",
    "-CABOND",
}

TEMPLATE_MODELS = {
    "SBM_AA": "-AA",
    "SBM_AA+gaussian": "-AAgaussian",
    "SBM_calpha": "-CA",
    "SBM_calpha+gaussian": "-CAgaussian",
    "SBM_2cg": "-AA2cg",
    "SBM_cr2": "-AAnbcr2",
    "SBM_AA+customContacts": "-AACC1",
    "SBM_AA+customContacts+customDihedrals": "-AACCD",
    "SBM_AA_BOND": "-AABOND",
    "SBM_AA_DIHE": "-AADIHE",
    "SBM_AA_DIHE4": "-AADIHE4",
    "SBM_match": "-AAMATCH",
}


def _has_model_flag(argv: list[str]) -> bool:
    return any(arg in MODEL_FLAGS for arg in argv)


def _template_flags(template: str, *, cg: bool = False) -> list[str]:
    path = Path(template.rstrip("/"))
    bif_path = _first_file(path, ".bif")
    nb_path = _first_file(path, ".nb")
    dynamic_template_flags: list[str] = []
    if bif_path is not None:
        dynamic_template_flags.extend(["-templateBif", str(bif_path)])
    if nb_path is not None:
        dynamic_template_flags.extend(["-templateNb", str(nb_path)])
    if cg:
        return ["-CA", *dynamic_template_flags, *_contact_flags_from_sif(path)]
    if path.name in TEMPLATE_MODELS:
        return [TEMPLATE_MODELS[path.name], *dynamic_template_flags]
    if path.name.startswith("SBM_CA"):
        return ["-CA", *dynamic_template_flags, *_contact_flags_from_sif(path)]
    if path.name.startswith("SBM_AA"):
        return ["-AA", *dynamic_template_flags, *_contact_flags_from_sif(path)]
    if path.name in {"temp.bifsif", "temp.cont.bifsif"}:
        model = "-AAMATCH" if (path / "CB.bif").exists() or (path / "comparelist").exists() else ("-CA" if "cont" in path.name else "-AA")
        flags = [model]
        flags.extend(_contact_flags_from_sif(path))
        if model == "-AAMATCH":
            flags.extend(_match_flags_from_nb(path))
        return flags
    return []


def _first_file(directory: Path, suffix: str) -> Path | None:
    if not directory.is_dir():
        return None
    return next(iter(sorted(directory.glob(f"*{suffix}"))), None)


def _contact_flags_from_sif(directory: Path) -> list[str]:
    sif = _first_file(directory, ".sif")
    if sif is None:
        return []
    try:
        root = ET.parse(sif).getroot()
    except ET.ParseError:
        return []
    contacts = root.find(".//Contacts")
    if contacts is None:
        return []
    method = contacts.attrib.get("method", "").strip()
    has_free_terms = any(
        (node.attrib.get("name", "").endswith("_free") or node.attrib.get("name") == "free")
        for node in root.findall(".//*")
    )
    has_gaussian_contacts = any(
        "gaussian" in node.attrib.get("name", "") or "gaussian" in node.attrib.get("func", "")
        for node in root.findall(".//*")
    )
    if method == "shadow" and has_free_terms:
        method = "shadow-free"
    if method == "cutoff" and has_gaussian_contacts:
        method = "cutoff-gaussian"
    if method not in {"shadow", "shadow-free", "cutoff", "cutoff-gaussian"}:
        method = "shadow" if method else ""
    flags: list[str] = []
    if method:
        flags.extend(["-contactMode", method])
    if contacts.attrib.get("contactDistance"):
        flags.extend(["-contactParam", contacts.attrib["contactDistance"]])
    if method in {"shadow", "shadow-free"}:
        if contacts.attrib.get("shadowRadius"):
            flags.extend(["-contactShadowSize", contacts.attrib["shadowRadius"]])
        if contacts.attrib.get("shadowRadiusBonded"):
            flags.extend(["-contactBondedRadius", contacts.attrib["shadowRadiusBonded"]])
    scaling = root.find(".//contactScaling[@name='stackingScale']") or root.find(".//contactScaling")
    if method in {"cutoff", "cutoff-gaussian"} and scaling is not None and scaling.attrib.get("scale"):
        flags.extend(["-contactStackScale", scaling.attrib["scale"]])
    dihedral_norm = root.find(".//dihedralNormalization")
    if dihedral_norm is not None and dihedral_norm.attrib.get("dihedralCounting"):
        flags.extend(["-dihedralCounting", dihedral_norm.attrib["dihedralCounting"]])
    return flags


def _match_flags_from_nb(directory: Path) -> list[str]:
    nb = _first_file(directory, ".nb")
    if nb is None:
        return []
    try:
        root = ET.parse(nb).getroot()
    except ET.ParseError:
        return []
    defaults = root.find(".//defaults")
    if defaults is None:
        return []
    flags: list[str] = []
    if defaults.attrib.get("gen-pairs"):
        flags.extend(["-matchGenPairs", defaults.attrib["gen-pairs"]])
    if defaults.attrib.get("fudgeLJ"):
        flags.extend(["-matchFudgeLJ", defaults.attrib["fudgeLJ"]])
    if defaults.attrib.get("fudgeQQ"):
        flags.extend(["-matchFudgeQQ", defaults.attrib["fudgeQQ"]])
    return flags


def translate_smogcheck_args(argv: list[str], *, stdin_text: str | None = None) -> list[str]:
    """Translate legacy SMOG-CHECK harness flags to native SMOG3 smog2 flags."""

    translated: list[str] = []
    template_flags: list[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "-opensmog":
            translated.append("-OpenSMOG")
            i += 1
            continue
        if arg in {"-SCMorig", "-keep4SCM"}:
            i += 1
            continue
        if arg in {"-t", "-tCG"}:
            if i + 1 >= len(argv):
                translated.append(arg)
                i += 1
                continue
            template = argv[i + 1]
            template_flags.extend(_template_flags(template, cg=(arg == "-tCG")))
            i += 2
            continue
        translated.append(arg)
        i += 1

    if any(arg in {"-v", "-help", "--help"} for arg in translated):
        return translated

    if not _has_model_flag(translated):
        template_model_flags = [flag for flag in template_flags if flag in MODEL_FLAGS]
        model_from_template = template_model_flags[-1] if template_model_flags else None
        has_gaussian_template = any(
            template_flags[idx] == "-contactMode" and idx + 1 < len(template_flags) and template_flags[idx + 1].endswith("-gaussian")
            for idx in range(len(template_flags))
        )
        if has_gaussian_template and model_from_template == "-AA":
            model_from_template = "-AAgaussian"
        elif has_gaussian_template and model_from_template == "-CA":
            model_from_template = "-CAgaussian"
        if model_from_template:
            translated.append(model_from_template)
        else:
            choice = stdin_text
            if choice is None and not sys.stdin.isatty():
                choice = sys.stdin.readline()
            choice = (choice or "").strip()
            if choice == "0":
                translated.append("-AA")
            elif choice == "2":
                translated.append("-CA")
    skip_value = False
    existing_option_names = {arg for arg in translated if arg.startswith("-")}
    for idx, flag in enumerate(template_flags):
        if skip_value:
            skip_value = False
            continue
        if flag in MODEL_FLAGS:
            continue
        if flag.startswith("-") and flag in existing_option_names and flag not in {
            "-templateBif",
            "-templateNb",
            "-contactMode",
            "-contactParam",
            "-contactShadowSize",
            "-contactBondedRadius",
            "-contactStackScale",
            "-dihedralCounting",
        }:
            skip_value = idx + 1 < len(template_flags) and not template_flags[idx + 1].startswith("-")
            continue
        translated.append(flag)
        if flag.startswith("-"):
            existing_option_names.add(flag)
    return translated


def _contact_output_path(argv: list[str]) -> Path | None:
    for idx, arg in enumerate(argv[:-1]):
        if arg == "-s":
            return Path(argv[idx + 1])
    return None


def _write_dropin_log(original: list[str], translated: list[str]) -> None:
    log_path = os.environ.get("SMOG3_DROPIN_LOG")
    if not log_path:
        return
    record = {
        "cwd": os.getcwd(),
        "pid": os.getpid(),
        "original": original,
        "translated": translated,
    }
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def _ensure_shadow_output(translated: list[str]) -> None:
    contacts = _contact_output_path(translated)
    if contacts is None or not contacts.exists():
        return
    shadow = contacts.with_suffix(contacts.suffix + ".ShadowOutput")
    if shadow.exists():
        return
    shutil.copyfile(contacts, shadow)


def main(argv: list[str] | None = None) -> int:
    original = list(sys.argv[1:] if argv is None else argv)
    os.environ["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    os.environ.setdefault("SMOG3_USE_SCM_DEFAULTS", "1")
    translated = translate_smogcheck_args(original)
    _write_dropin_log(original, translated)
    rc = smog2_native_main(translated)
    if rc == 0:
        _ensure_shadow_output(translated)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
