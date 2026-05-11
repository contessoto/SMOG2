from pathlib import Path

from smog3 import smogcheck_parity


def test_parse_smogcheck_testlist_discovers_numbered_cases():
    cases = smogcheck_parity.parse_testlist()
    assert len(cases) == 115
    assert cases[0].case_id == 1
    assert cases[-1].case_id == 115
    assert sum(1 for case in cases if case.opensmog) == 14


def test_case_translation_for_opensmog_and_user_contacts(tmp_path: Path):
    cases = {case.case_id: case for case in smogcheck_parity.parse_testlist()}
    out = smogcheck_parity.ROOT / "parity_runs" / "unit" / "case94"
    args94 = smogcheck_parity._candidate_args(cases[94], out)
    assert args94 is not None
    assert "-OpenSMOG" in args94
    assert str(out / "model.xml") not in args94

    args50 = smogcheck_parity._candidate_args(cases[50], out)
    assert args50 is not None
    assert "-c" in args50
    user_contact_index = len(args50) - 1 - list(reversed(args50)).index("-c")
    assert "2ci2_v2.contacts" in args50[user_contact_index + 1]


def test_unsupported_cases_are_classified_by_feature_group():
    cases = {case.case_id: case for case in smogcheck_parity.parse_testlist()}
    assert smogcheck_parity._supported(cases[41])[0] is True
    assert smogcheck_parity.feature_group(cases[41]) == "CA coarse-graining"
    assert smogcheck_parity._supported(cases[53])[0] is True
    assert smogcheck_parity.feature_group(cases[53]) == "template/map variants"
    assert smogcheck_parity._supported(cases[114])[0] is True
    assert smogcheck_parity.feature_group(cases[114]) == "freecoor/interactive"


def test_interactive_cases_use_reproducible_equivalent_model_flags():
    cases = {case.case_id: case for case in smogcheck_parity.parse_testlist()}
    out = smogcheck_parity.ROOT / "parity_runs" / "unit" / "interactive"

    assert smogcheck_parity._baseline_model_flags(cases[114]) == ["-AA"]
    assert smogcheck_parity._candidate_args(cases[114], out)[-1] == "-AA"  # type: ignore[index]
    assert smogcheck_parity._baseline_model_flags(cases[115]) == ["-CA"]
    assert smogcheck_parity._candidate_args(cases[115], out)[-1] == "-CA"  # type: ignore[index]
