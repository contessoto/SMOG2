from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .parity_direct import compare_existing_dirs
ROOT = Path(__file__).resolve().parents[2]
TESTLIST = ROOT / "SMOG-CHECK" / "share" / "settings" / "smog2.testlist"
DOCKER_SMOGCHECK_TEMPLATES = "/workdir/SMOG-CHECK/share/templates"

KNOWN_MODELS = {
    "AA",
    "CA",
    "AA-2cg",
    "AA-nb-cr2",
    "AA-match",
    "AA-CC1",
    "AA-CCD",
    "AA-BOND",
    "AA-DIHE",
    "AA-DIHE4",
    "AA-interactive",
    "CA-interactive",
}

MODEL_FLAGS = {
    ("AA", False): ["-AA"],
    ("AA", True): ["-AAgaussian"],
    ("CA", False): ["-CA"],
    ("CA", True): ["-CAgaussian"],
    ("AA-2cg", False): ["-AA2cg"],
    ("AA-nb-cr2", False): ["-AAnbcr2"],
    ("AA-CC1", False): ["-AACC1"],
    ("AA-CCD", False): ["-AACCD"],
    ("AA-BOND", False): ["-AABOND"],
    ("AA-DIHE", False): ["-AADIHE"],
    ("AA-DIHE4", False): ["-AADIHE4"],
    ("AA-match", False): ["-AAMATCH"],
    ("AA-interactive", False): ["-AA"],
    ("CA-interactive", False): ["-CA"],
}

BASELINE_TEMPLATE_FLAGS = {
    "AA-2cg": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_2cg"],
    "AA-nb-cr2": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_cr2"],
    "AA-CC1": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA+customContacts"],
    "AA-CCD": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA+customContacts+customDihedrals"],
    "AA-BOND": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA_BOND"],
    "AA-DIHE": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA_DIHE"],
    "AA-DIHE4": [f"-t", f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA_DIHE4"],
}


@dataclass(frozen=True)
class SmogcheckCase:
    case_id: int
    stem: str
    modifiers: tuple[str, ...]
    model: str
    contact_model: str
    params: tuple[str, ...]
    raw: str

    @property
    def pdb(self) -> str:
        return f"SMOG-CHECK/share/PDB.files/{self.stem}.pdb"

    @property
    def opensmog(self) -> bool:
        return "OpenSMOG" in self.modifiers

    @property
    def freecoor(self) -> bool:
        return "freecoor" in self.modifiers

    @property
    def interactive(self) -> bool:
        return self.model.endswith("-interactive")

    @property
    def gaussian(self) -> bool:
        return "gaussian" in self.contact_model

    @property
    def user_contacts(self) -> bool:
        return "userC" in self.contact_model


def parse_testlist(path: Path = TESTLIST) -> list[SmogcheckCase]:
    cases: list[SmogcheckCase] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";") or ";" not in line:
            continue
        left, case_txt = line.rsplit(";", 1)
        try:
            case_id = int(case_txt.strip())
        except ValueError:
            continue
        tokens = left.split()
        if len(tokens) < 3:
            continue
        stem = tokens[0]
        model_index = next((idx for idx, token in enumerate(tokens[1:], start=1) if token in KNOWN_MODELS), None)
        if model_index is None or model_index + 1 >= len(tokens):
            continue
        cases.append(
            SmogcheckCase(
                case_id=case_id,
                stem=stem,
                modifiers=tuple(tokens[1:model_index]),
                model=tokens[model_index],
                contact_model=tokens[model_index + 1],
                params=tuple(tokens[model_index + 2 :]),
                raw=line,
            )
        )
    return cases


def feature_group(case: SmogcheckCase) -> str:
    if case.interactive or case.freecoor:
        return "freecoor/interactive"
    if case.opensmog:
        return "OpenSMOG XML"
    if case.model == "CA" or case.model == "CA-interactive":
        return "CA coarse-graining"
    if case.contact_model == "shadow-free":
        return "shadow-free/custom topology"
    if case.contact_model in {"shadow", "cutoff", "cutoff-gaussian", "shadow-match"}:
        return "template/map variants"
    if case.user_contacts:
        return "user contacts"
    if case.model not in {"AA", "CA"}:
        return "template/map variants"
    return "topology parameter parity"


def _rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _shell_join(args: list[str]) -> str:
    return " ".join(subprocess.list2cmdline([arg]) for arg in args)


def _model_flags(case: SmogcheckCase) -> list[str] | None:
    return MODEL_FLAGS.get((case.model, case.gaussian))


def _baseline_model_flags(case: SmogcheckCase) -> list[str] | None:
    if case.interactive:
        return _model_flags(case)
    if case.contact_model in {"shadow", "shadow-free", "cutoff", "cutoff-gaussian"}:
        if case.model == "AA":
            return ["-t", "temp.bifsif/"]
        if case.model == "CA":
            return ["-tCG", "temp.bifsif/", "-t", "temp.cont.bifsif"]
    if case.contact_model == "shadow-match" and case.model == "AA-match":
        return ["-t", "temp.bifsif/"]
    if case.contact_model == "default" and case.model in BASELINE_TEMPLATE_FLAGS:
        return list(BASELINE_TEMPLATE_FLAGS[case.model])
    return _model_flags(case)


def _param(case: SmogcheckCase, index: int, default: str) -> str:
    return case.params[index] if index < len(case.params) else default


def _fmt_num(value: float) -> str:
    return f"{value:.12g}"


def _nondefault_values(case: SmogcheckCase) -> dict[str, str]:
    params = case.params
    if case.contact_model in {"shadow", "shadow-free"}:
        return {
            "CONTD": _param(case, 0, "5.0"),
            "CONTR": _param(case, 1, "1.0"),
            "STACKSCALE": "1.0",
            "R_CD": _param(case, 2, "1.2"),
            "R_P_BB_SC": _param(case, 3, "1.0"),
            "R_N_SC_BB": _param(case, 4, "2.0"),
            "PRO_DIH": _param(case, 5, "1.0"),
            "NA_DIH": _param(case, 6, "1.0"),
            "LIGAND_DIH": _param(case, 7, "1.0"),
            "sigma": _param(case, 8, "2.5"),
            "epsilon": _param(case, 9, "0.01"),
            "epsilonCAC": _param(case, 10, "1.0"),
            "epsilonCAD": _param(case, 11, "1.4"),
            "sigmaCA": _param(case, 12, "5.0"),
            "massNB2": _param(case, 13, "0.2"),
            "chargeNB2": _param(case, 14, "-1"),
            "C6NB2": _param(case, 15, "1E-6"),
            "C12NB2": _param(case, 16, "3E-9"),
            "chargeAT": _param(case, 17, "1.0"),
            "DIHEDCOUNT": "0",
        }
    if case.contact_model in {"cutoff", "cutoff-gaussian"}:
        return {
            "CONTD": _param(case, 0, "5.0"),
            "CONTR": "0.0",
            "STACKSCALE": _param(case, 1, "1.0"),
            "R_CD": _param(case, 2, "1.2"),
            "R_P_BB_SC": _param(case, 3, "1.0"),
            "R_N_SC_BB": _param(case, 4, "2.0"),
            "PRO_DIH": _param(case, 5, "1.0"),
            "NA_DIH": _param(case, 6, "1.0"),
            "LIGAND_DIH": _param(case, 7, "1.0"),
            "sigma": _param(case, 8, "2.5"),
            "epsilon": _param(case, 9, "0.01"),
            "epsilonCAC": _param(case, 10, "1.0"),
            "epsilonCAD": _param(case, 11, "1.4"),
            "sigmaCA": _param(case, 12, "5.0"),
            "massNB2": _param(case, 13, "0.2"),
            "chargeNB2": _param(case, 14, "-1"),
            "C6NB2": _param(case, 15, "1E-6"),
            "C12NB2": _param(case, 16, "3E-9"),
            "chargeAT": _param(case, 17, "1.0"),
            "DIHEDCOUNT": _param(case, 18, "0"),
        }
    return {}


def _aa_template_setup(case: SmogcheckCase) -> str:
    values = _nondefault_values(case)
    if not values:
        return ""
    template = f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA"
    pro_dih = float(values["PRO_DIH"])
    na_dih = float(values["NA_DIH"])
    parm_p_bb = pro_dih
    parm_p_sc = pro_dih / float(values["R_P_BB_SC"])
    parm_n_bb = na_dih
    parm_n_sc = na_dih * float(values["R_N_SC_BB"])
    sigma_nm = float(values["sigma"]) / 10.0
    c12_nb1 = (sigma_nm**12) * float(values["epsilon"])
    sif_name = {
        "shadow": "AA-test.shadow.sif",
        "shadow-free": "AA-test.shadow.free.sif",
        "cutoff": "AA-test.cutoff.sif",
        "cutoff-gaussian": "AA-test.cutoff.gaussian.sif",
    }[case.contact_model]
    if case.contact_model == "shadow-free":
        nb_name = "AA-test.free.nb"
        bif_name = "AA-test.free.bif"
        b_name = "AA-test.free.b"
    else:
        nb_name = "AA-test.gaussian.nb" if case.contact_model == "cutoff-gaussian" else "AA-test.nb"
        bif_name = "AA-test.bif"
        b_name = "AA-test.b"
    sif_replacements = (
        f"s/PARM_C_D/{values['R_CD']}/g;"
        f"s/PARM_P_BB/{_fmt_num(parm_p_bb)}/g;"
        f"s/PARM_P_SC/{_fmt_num(parm_p_sc)}/g;"
        f"s/PARM_N_BB/{_fmt_num(parm_n_bb)}/g;"
        f"s/PARM_N_SC/{_fmt_num(parm_n_sc)}/g;"
        f"s/CUTDIST/{values['CONTD']}/g;"
    )
    if case.contact_model in {"shadow", "shadow-free"}:
        sif_replacements += f"s/SCM_R/{values['CONTR']}/g;s/SCM_BR/0.5/g;"
    else:
        sif_replacements += f"s/RESCALE/{values['STACKSCALE']}/g;s/DIHEDCOUNT/{values['DIHEDCOUNT']}/g;"
    sif_replacements += "s/MINVERSION/2.4.5/g"
    nb_replacements = (
        f"s/PARM_MASS/{values['massNB2']}/g;"
        f"s/PARM_chargeNB/{values['chargeNB2']}/g;"
        f"s/PARM_C6_2/{values['C6NB2']}/g;"
        f"s/PARM_C12_2/{values['C12NB2']}/g;"
        f"s/PARM_C12/{_fmt_num(c12_nb1)}/g"
    )
    return f"""
mkdir -p temp.bifsif
sed "{sif_replacements}" {template}/{sif_name} > temp.bifsif/tmp.sif
sed "{nb_replacements}" {template}/{nb_name} > temp.bifsif/tmp.nb
cp {template}/{bif_name} temp.bifsif/tmp.bif
cp {template}/{b_name} temp.bifsif/tmp.b
cp {template}/extras temp.bifsif/test.extras
"""


def _ca_template_setup(case: SmogcheckCase) -> str:
    values = _nondefault_values(case)
    if not values:
        return ""
    ca_template = f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_calpha"
    aa_template = f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_AA_STATIC"
    sigma_ca_nm = float(values["sigmaCA"]) / 10.0
    c12_ca = sigma_ca_nm**12
    epsilon_cad3 = float(values["epsilonCAD"]) / 2.0
    if case.contact_model == "shadow":
        contact_sif = (
            f'sed "s/CUTDIST/{values["CONTD"]}/g;s/SCM_R/{values["CONTR"]}/g;'
            f's/SCM_BR/0.5/g;s/MINVERSION/2.4.5/g" '
            f"{aa_template}/AA-test.shadow.sif > temp.cont.bifsif/tmp.cont.sif"
        )
    else:
        contact_sif = (
            f'sed "s/CUTDIST/{values["CONTD"]}/g;s/RESCALE/{values["STACKSCALE"]}/g;'
            f's/DIHEDCOUNT/{values["DIHEDCOUNT"]}/g;s/MINVERSION/2.4.5/g" '
            f"{aa_template}/AA-test.cutoff.sif > temp.cont.bifsif/tmp.cont.sif"
        )
    return f"""
mkdir -p temp.bifsif temp.cont.bifsif
sed "s/EPS_CONT/{values['epsilonCAC']}/g;s/EPS_DIH/{values['epsilonCAD']}/g;s/EPS_dih3/{_fmt_num(epsilon_cad3)}/g;s/MINVERSION/2.4.5/g" {ca_template}/CA-test.sif > temp.bifsif/tmp.sif
sed "s/PARM_C12/{_fmt_num(c12_ca)}/g;s/EPS_CONT/{values['epsilonCAC']}/g" {ca_template}/CA-test.nb > temp.bifsif/tmp.nb
sed "s/EPS_CONT/{values['epsilonCAC']}/g;s/EPS_DIH/{values['epsilonCAD']}/g;s/EPS_dih3/{_fmt_num(epsilon_cad3)}/g" {ca_template}/CA-test.b > temp.bifsif/tmp.b
cp {ca_template}/CA-test.bif temp.bifsif/tmp.bif
cp {aa_template}/AA-test.bif temp.cont.bifsif/tmp.cont.bif
cp {aa_template}/AA-test.nb temp.cont.bifsif/tmp.cont.nb
cp {aa_template}/AA-test.b temp.cont.bifsif/tmp.cont.b
{contact_sif}
"""


def _match_template_setup(case: SmogcheckCase) -> str:
    genpairs = _param(case, 2, "0")
    fudgelj = _param(case, 3, "1.0")
    fudgeqq = _param(case, 4, "1.0")
    template = f"{DOCKER_SMOGCHECK_TEMPLATES}/SBM_match"
    return f"""
mkdir -p temp.bifsif
cp -r {template}/. temp.bifsif/
sed "s/GENPAIRS/{genpairs}/g;s/FUDGELJ/{fudgelj}/g;s/FUDGEQQ/{fudgeqq}/g" {template}/CB.nb > temp.bifsif/CB.nb
"""


def _supported(case: SmogcheckCase) -> tuple[bool, str]:
    if _model_flags(case) is None:
        return False, f"model {case.model} is not implemented natively"
    return True, ""


def _baseline_setup(case: SmogcheckCase) -> str:
    if case.contact_model in {"shadow", "shadow-free", "cutoff", "cutoff-gaussian"}:
        if case.model == "AA":
            return _aa_template_setup(case)
        if case.model == "CA":
            return _ca_template_setup(case)
    if case.contact_model == "shadow-match" and case.model == "AA-match":
        return _match_template_setup(case)
    return ""


def _baseline_args(case: SmogcheckCase, outdir: Path) -> list[str] | None:
    flags = _baseline_model_flags(case)
    if flags is None:
        return None
    args = [
        "smog2",
        "-i",
        f"/workdir/{case.pdb}",
        "-SCMorig",
        "-keep4SCM",
        "-o",
        f"/workdir/{_rel(outdir / 'model.top')}",
        "-g",
        f"/workdir/{_rel(outdir / 'model.gro')}",
        "-n",
        f"/workdir/{_rel(outdir / 'model.ndx')}",
        "-s",
        f"/workdir/{_rel(outdir / 'model.contacts')}",
    ]
    if case.opensmog:
        args.extend(["-OpenSMOG", "-OpenSMOGxml", f"/workdir/{_rel(outdir / 'model.xml')}"])
    if case.freecoor:
        args.append("-freecoor")
    args.extend(flags)
    if case.user_contacts:
        args.extend(["-c", f"/workdir/SMOG-CHECK/share/PDB.files/{case.stem}.contacts"])
    return args


def _candidate_args(case: SmogcheckCase, outdir: Path) -> list[str] | None:
    flags = _model_flags(case)
    if flags is None:
        return None
    args = [
        "python3",
        "-c",
        "from smog3.smog2_native import main; import sys; raise SystemExit(main(sys.argv[1:]))",
        "-i",
        case.pdb,
        "-o",
        _rel(outdir / "model.top"),
        "-g",
        _rel(outdir / "model.gro"),
        "-n",
        _rel(outdir / "model.ndx"),
        "-s",
        _rel(outdir / "model.contacts"),
    ]
    if case.opensmog:
        args.extend(["-OpenSMOG", "-OpenSMOGxml", _rel(outdir / "model.xml")])
    if case.freecoor:
        args.append("-freecoor")
    args.extend(flags)
    if case.user_contacts:
        args.extend(["-c", f"SMOG-CHECK/share/PDB.files/{case.stem}.contacts"])
    if case.model == "AA-match":
        args.extend(
            [
                "-matchGenPairs",
                _param(case, 2, "0"),
                "-matchFudgeLJ",
                _param(case, 3, "1.0"),
                "-matchFudgeQQ",
                _param(case, 4, "1.0"),
            ]
        )
    if case.contact_model in {"shadow", "shadow-free", "shadow-match"}:
        native_mode = "shadow-free" if case.contact_model == "shadow-free" else "shadow"
        args.extend(
            [
                "-contactMode",
                native_mode,
                "-contactParam",
                _param(case, 0, "5.0"),
                "-contactShadowSize",
                _param(case, 1, "1.0"),
                "-contactBondedRadius",
                "0.5",
            ]
        )
    elif case.contact_model in {"cutoff", "cutoff-gaussian"}:
        args.extend(
            [
                "-contactMode",
                case.contact_model,
                "-contactParam",
                _param(case, 0, "5.0"),
                "-contactStackScale",
                _param(case, 1, "1.0"),
                "-dihedralCounting",
                _param(case, 18, "0"),
            ]
        )
    return args


def _run_baseline(case: SmogcheckCase, outdir: Path, image: str) -> tuple[int, str, list[str]]:
    args = _baseline_args(case, outdir)
    if args is None:
        return 99, "unsupported baseline argument mapping", []
    setup = _baseline_setup(case)
    command = _shell_join(args)
    if case.interactive:
        choice = "2" if case.model == "CA-interactive" else "0"
        command = f"printf '%s\\n' {subprocess.list2cmdline([choice])} | {command}"
    script = f"""
set -euo pipefail
cd /workdir
cd {_rel(outdir)}
{setup}
    {command}
"""
    proc = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{ROOT}:/workdir", image, "bash", "-lc", script],
        capture_output=True,
        text=True,
    )
    return proc.returncode, proc.stdout + proc.stderr, args


def _run_candidate(case: SmogcheckCase, outdir: Path) -> tuple[int, str, list[str]]:
    args = _candidate_args(case, outdir)
    if args is None:
        return 99, "unsupported candidate argument mapping", []
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    env["SMOG3_LEGACY_PERL_FALLBACK"] = "0"
    env["SMOG3_USE_SCM_DEFAULTS"] = "1"
    proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout + proc.stderr, args


def _line_counts(directory: Path, include_xml: bool) -> dict[str, int | None]:
    names = ["model.top", "model.gro", "model.ndx", "model.contacts"] + (["model.xml"] if include_xml else [])
    return {name: (len((directory / name).read_text().splitlines()) if (directory / name).exists() else None) for name in names}


def _case_status(baseline_rc: int, candidate_rc: int, comparison_ok: bool) -> str:
    if baseline_rc != 0:
        return "BASELINE_ERROR"
    if candidate_rc != 0:
        return "CANDIDATE_ERROR"
    return "PASS" if comparison_ok else "FAIL"


def _select_cases(all_cases: list[SmogcheckCase], ns: argparse.Namespace) -> list[SmogcheckCase]:
    if ns.cases:
        wanted = {int(x) for x in ns.cases.split(",") if x.strip()}
        selected = [case for case in all_cases if case.case_id in wanted]
    elif ns.range:
        lo_txt, hi_txt = ns.range.split("-", 1)
        lo, hi = int(lo_txt), int(hi_txt)
        selected = [case for case in all_cases if lo <= case.case_id <= hi]
    elif ns.all:
        selected = list(all_cases)
    else:
        raise SystemExit("Choose --cases, --range, or --all")
    if ns.stop_after is not None:
        selected = selected[: ns.stop_after]
    return selected


def run_campaign(cases: list[SmogcheckCase], out_root: Path, image: str) -> dict:
    if shutil.which("docker") is None:
        return {"ok": False, "reason": "docker not available", "cases": []}
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True)
    reports_dir = out_root / "reports"
    reports_dir.mkdir()

    report: dict = {
        "ok": True,
        "image": image,
        "total_selected": len(cases),
        "cases": [],
        "summary": {},
        "feature_groups": {},
    }
    for case in cases:
        supported, reason = _supported(case)
        case_root = out_root / f"case{case.case_id}"
        baseline_dir = case_root / "baseline"
        candidate_dir = case_root / "candidate"
        baseline_dir.mkdir(parents=True)
        candidate_dir.mkdir(parents=True)
        group = feature_group(case)
        if not supported:
            entry = {
                "case": case.case_id,
                "pdb": case.stem,
                "model": case.model,
                "contact_model": case.contact_model,
                "feature_group": group,
                "status": "UNSUPPORTED",
                "reason": reason,
                "raw": case.raw,
            }
            report["cases"].append(entry)
            (reports_dir / f"case{case.case_id}.json").write_text(json.dumps(entry, indent=2))
            report["ok"] = False
            print(f"case {case.case_id}: UNSUPPORTED ({group})", flush=True)
            continue

        brc, bout, bargs = _run_baseline(case, baseline_dir, image)
        crc, cout, cargs = _run_candidate(case, candidate_dir)
        comparison = compare_existing_dirs(baseline_dir, candidate_dir, include_xml=case.opensmog)
        status = _case_status(brc, crc, comparison["ok"])
        entry = {
            "case": case.case_id,
            "pdb": case.stem,
            "model": case.model,
            "contact_model": case.contact_model,
            "feature_group": group,
            "status": status,
            "baseline_rc": brc,
            "candidate_rc": crc,
            "baseline_args": bargs,
            "candidate_args": cargs,
            "baseline_tail": bout[-2000:],
            "candidate_tail": cout[-2000:],
            "baseline_line_counts": _line_counts(baseline_dir, case.opensmog),
            "candidate_line_counts": _line_counts(candidate_dir, case.opensmog),
            "comparisons": comparison["comparisons"],
            "raw": case.raw,
        }
        report["cases"].append(entry)
        (reports_dir / f"case{case.case_id}.json").write_text(json.dumps(entry, indent=2))
        report["ok"] = bool(report["ok"] and status == "PASS")
        print(f"case {case.case_id}: {status} ({group})", flush=True)

    summary: dict[str, int] = {}
    groups: dict[str, dict[str, int]] = {}
    for entry in report["cases"]:
        status = entry["status"]
        summary[status] = summary.get(status, 0) + 1
        group = entry["feature_group"]
        groups.setdefault(group, {})
        groups[group][status] = groups[group].get(status, 0) + 1
    report["summary"] = summary
    report["feature_groups"] = groups
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--cases")
    group.add_argument("--range")
    group.add_argument("--all", action="store_true")
    parser.add_argument("--stop-after", type=int, default=None)
    parser.add_argument("--out-root", default="parity_runs/all")
    parser.add_argument("--report-json", default="parity_all_summary.json")
    parser.add_argument("--image", default="smogserver/smog2:stable")
    ns = parser.parse_args(argv)

    all_cases = parse_testlist()
    selected = _select_cases(all_cases, ns)
    report = run_campaign(selected, ROOT / ns.out_root, ns.image)
    report["discovered"] = {
        "total_cases": len(all_cases),
        "direct_smog2_cases": len(all_cases),
        "helper_tool_cases": 0,
    }
    (ROOT / ns.report_json).write_text(json.dumps(report, indent=2))

    print("case status feature")
    for entry in report.get("cases", []):
        print(f"{entry['case']:>4} {entry['status']:<15} {entry['feature_group']}")
    print("\nsummary")
    for status in ["PASS", "FAIL", "UNSUPPORTED", "BASELINE_ERROR", "CANDIDATE_ERROR"]:
        print(f"  {status:<15} {report.get('summary', {}).get(status, 0)}")
    print(f"\nreport written to {ns.report_json}")
    return 0 if report.get("ok") else 2


if __name__ == "__main__":
    raise SystemExit(main())
