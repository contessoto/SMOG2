from pathlib import Path
from types import SimpleNamespace

from smog3 import selected_parity, smog2_native


def test_selected_case_definitions_cover_requested_cases():
    assert set(selected_parity.SELECTED_CASES) == {1, 21, 41, 50, 56, 94}
    assert selected_parity.SELECTED_CASES[50].baseline_args[:2] == ("-AA", "-c")
    assert selected_parity.SELECTED_CASES[50].candidate_args[:2] == ("-AA", "-c")
    assert selected_parity.SELECTED_CASES[56].baseline_setup
    assert selected_parity.SELECTED_CASES[94].include_xml is True


def test_selected_parity_reports_per_case_comparisons(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(selected_parity.shutil, "which", lambda name: "/usr/bin/docker" if name == "docker" else None)

    def fake_baseline(case, outdir, image):
        for name in ["model.top", "model.gro", "model.ndx", "model.contacts"]:
            (outdir / name).write_text("same\n")
        if case.include_xml:
            (outdir / "model.xml").write_text("same\n")
        return 0, "baseline ok", ["smog2"]

    def fake_candidate(case, outdir):
        for name in ["model.top", "model.gro", "model.ndx", "model.contacts"]:
            (outdir / name).write_text("same\n")
        if case.include_xml:
            (outdir / "model.xml").write_text("same\n")
        return 0, "candidate ok", ["python3"]

    monkeypatch.setattr(selected_parity, "_run_baseline", fake_baseline)
    monkeypatch.setattr(selected_parity, "_run_candidate", fake_candidate)

    report = selected_parity.run_selected([1, 94], tmp_path / "selected", "smogserver/smog2:stable")
    assert report["ok"] is True
    assert [case["case"] for case in report["cases"]] == [1, 94]
    assert report["cases"][1]["comparisons"]["model.xml"]["match"] is True


def test_native_accepts_smog2_c_alias_for_user_contacts(tmp_path: Path):
    pdb = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files" / "2ci2_v2.pdb"
    user_contacts = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files" / "2ci2_v2.contacts"
    base = tmp_path / "case50"
    rc = smog2_native.main([
        "-i", str(pdb),
        "-AA",
        "-c", str(user_contacts),
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])
    assert rc == 0
    assert base.with_suffix(".contacts").read_text() == user_contacts.read_text()


def test_selected_aa_gaussian_uses_full_scm_topology(monkeypatch, tmp_path: Path):
    pdb = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files" / "1A01-AMP.pdb"
    base = tmp_path / "case21"

    def fake_scm(_coord, _top, out_contacts, _atoms):
        out_contacts.write_text("1 1 2 1195\n")
        return [(1, 1, ("2", "1195"))]

    monkeypatch.setenv("SMOG3_USE_SCM_DEFAULTS", "1")
    monkeypatch.setattr(smog2_native, "_generate_contacts_with_scm", fake_scm)

    rc = smog2_native.main([
        "-i", str(pdb),
        "-AAgaussian",
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])

    top_text = base.with_suffix(".top").read_text()
    assert rc == 0
    assert "[ angles ]" in top_text
    assert "[ dihedrals ]" in top_text
    assert "\t6\t" in top_text


def test_selected_shadow_free_uses_java_scm_parameters(monkeypatch, tmp_path: Path):
    pdb = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files" / "1AKEapo_v2.pdb"
    base = tmp_path / "case56"
    seen = {}

    def fake_scm(_coord, _top, out_contacts, _atoms, **kwargs):
        seen.update(kwargs)
        out_contacts.write_text("1 1 3 584\n")
        return [(1, 1, ("3", "584"))]

    monkeypatch.setenv("SMOG3_USE_SCM_DEFAULTS", "1")
    monkeypatch.setattr(smog2_native, "_generate_contacts_with_scm", fake_scm)

    rc = smog2_native.main([
        "-i", str(pdb),
        "-AA",
        "-contactMode", "shadow-free",
        "-contactParam", "5.0",
        "-o", str(base.with_suffix(".top")),
        "-g", str(base.with_suffix(".gro")),
        "-n", str(base.with_suffix(".ndx")),
        "-s", str(base.with_suffix(".contacts")),
    ])

    assert rc == 0
    assert seen["mode"] == "shadow"
    assert seen["cutoff"] == 5.0
    assert seen["shadow_size"] == 1.4
