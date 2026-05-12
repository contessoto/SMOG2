from __future__ import annotations

from pathlib import Path

from smog3 import tutorial_validation


def test_tutorial_manifest_has_representative_implemented_cases() -> None:
    cases = {case.case_id: case for case in tutorial_validation.tutorial_cases()}

    for case_id in {
        "standard_aa_ci2",
        "standard_ca_ci2",
        "opensmog_aa_ci2",
        "opensmog_ca_ci2",
        "aa_1a01_amp_ligand",
        "aa_dna_terminal",
        "aa_trna",
        "aa_disulfide_local",
        "ca_disulfide_local",
        "aa_user_contacts_ci2",
    }:
        assert cases[case_id].implemented

    assert any(case.status == "manual_input_required" for case in cases.values())
    assert any(case.status == "not_generation_test" for case in cases.values())
    assert "MISSING_DOWNLOAD" in tutorial_validation.ALL_STATUSES


def test_tutorial_opensmog_cases_use_supported_xml_flags() -> None:
    cases = {case.case_id: case for case in tutorial_validation.tutorial_cases()}
    aa = cases["opensmog_aa_ci2"]
    ca = cases["opensmog_ca_ci2"]

    assert aa.include_xml
    assert ca.include_xml
    assert aa.raw["smog3_args"] == ["-AA", "-OpenSMOG", "-OpenSMOGxml", "{model_xml}"]
    assert ca.raw["smog3_args"] == ["-CA", "-OpenSMOG", "-OpenSMOGxml", "{model_xml}"]


def test_tutorial_command_builder_separates_docker_and_local_paths(tmp_path: Path) -> None:
    case = {case.case_id: case for case in tutorial_validation.tutorial_cases()}["aa_user_contacts_ci2"]
    context = {
        "input_pdb_rel": "validation/tutorials/runs/test/inputs/aa_user_contacts_ci2/2ci2_v2.pdb",
        "input_pdb_docker": "/workdir/validation/tutorials/runs/test/inputs/aa_user_contacts_ci2/2ci2_v2.pdb",
        "user_contacts_rel": "validation/tutorials/runs/test/inputs/aa_user_contacts_ci2/2ci2_v2.contacts",
        "user_contacts_docker": "/workdir/validation/tutorials/runs/test/inputs/aa_user_contacts_ci2/2ci2_v2.contacts",
    }
    baseline_dir = tutorial_validation.ROOT / "validation/tutorials/runs/test/smog2_baseline/aa_user_contacts_ci2"
    candidate_dir = tutorial_validation.ROOT / "validation/tutorials/runs/test/smog3_candidate/aa_user_contacts_ci2"

    smog2_args = tutorial_validation._base_smog2_args(case, context, baseline_dir)
    smog3_args = tutorial_validation._base_smog3_args(case, context, candidate_dir, use_installed=False)

    assert "/workdir/validation/tutorials/runs/test/inputs/aa_user_contacts_ci2/2ci2_v2.contacts" in smog2_args
    assert "validation/tutorials/runs/test/inputs/aa_user_contacts_ci2/2ci2_v2.contacts" in smog3_args
    assert smog3_args[:3] == ["python3", "-m", "smog3.smog2_native"]


def test_tutorial_no_perl_sentinel_records_attempt(tmp_path: Path, monkeypatch) -> None:
    no_perl_bin = tmp_path / "bin"
    perl_log = tmp_path / "perl.log"
    tutorial_validation._write_no_perl_sentinel(no_perl_bin)
    monkeypatch.setenv("SMOG3_PERL_SENTINEL_LOG", str(perl_log))

    import subprocess

    result = subprocess.run([str(no_perl_bin / "perl"), "-e", "print 1"], capture_output=True, text=True)

    assert result.returncode == 127
    assert "SMOG3 attempted to invoke perl" in perl_log.read_text()
