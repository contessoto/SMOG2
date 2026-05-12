from pathlib import Path
import re

from smog3 import smog2_native
from smog3 import smogcheck_dropin_smog2 as dropin


ROOT = Path(__file__).resolve().parents[1]


def _top_section_rows(text: str, section: str) -> list[list[str]]:
    rows: list[list[str]] = []
    in_section = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == f"[ {section} ]":
            in_section = True
            continue
        if in_section and stripped.startswith("["):
            break
        if in_section and stripped and not stripped.startswith(";"):
            rows.append(stripped.split())
    return rows


def test_dropin_translates_template_and_scm_harness_flags() -> None:
    args = dropin.translate_smogcheck_args(
        [
            "-SCMorig",
            "-keep4SCM",
            "-i",
            "share/PDB.files/1A01-AMP.pdb",
            "-t",
            "/tmp/smog2-layout/SBM_AA",
            "-s",
            "model.contacts",
        ]
    )

    assert "-SCMorig" not in args
    assert "-keep4SCM" not in args
    assert "-t" not in args
    assert args[-1] == "-AA"


def test_dropin_translates_cg_template_to_ca_model() -> None:
    args = dropin.translate_smogcheck_args(
        [
            "-i",
            "share/PDB.files/1A01-AMP.pdb",
            "-tCG",
            "/tmp/smog2-layout/SBM_calpha",
        ]
    )

    assert args[-1] == "-CA"


def test_dropin_translates_generated_match_template(tmp_path: Path) -> None:
    template = tmp_path / "temp.bifsif"
    template.mkdir()
    (template / "CB.bif").write_text("<bif />\n", encoding="utf-8")
    (template / "CB.sif").write_text(
        """<?xml version='1.0'?>
<sif><settings>
<Groups><groupRatios contacts="1" dihedrals="1"/></Groups>
<Contacts method="shadow" contactDistance="10" shadowRadius="2" shadowRadiusBonded="0.5"/>
<dihedralNormalization dihedralCounting="0"/>
</settings></sif>
""",
        encoding="utf-8",
    )
    (template / "CB.nb").write_text(
        """<?xml version='1.0'?>
<nb><defaults gen-pairs="1" nbfunc="1" gmx-combination-rule="1" fudgeLJ="1.2" fudgeQQ="1.5"/></nb>
""",
        encoding="utf-8",
    )

    args = dropin.translate_smogcheck_args(["-i", "x.pdb", "-t", str(template)])

    assert "-AAMATCH" in args
    assert args[args.index("-contactMode") + 1] == "shadow"
    assert args[args.index("-contactParam") + 1] == "10"
    assert args[args.index("-matchGenPairs") + 1] == "1"
    assert args[args.index("-matchFudgeLJ") + 1] == "1.2"
    assert args[args.index("-matchFudgeQQ") + 1] == "1.5"


def test_dropin_detects_shadow_free_template(tmp_path: Path) -> None:
    template = tmp_path / "temp.bifsif"
    template.mkdir()
    (template / "tmp.sif").write_text(
        """<?xml version='1.0'?>
<sif>
<functions><function name="angle_free" directive="angles"/></functions>
<settings><Contacts method="shadow" contactDistance="5.0" shadowRadius="1.4" shadowRadiusBonded="0.5"/></settings>
</sif>
""",
        encoding="utf-8",
    )

    args = dropin.translate_smogcheck_args(["-i", "x.pdb", "-t", str(template)])

    assert "-AA" in args
    assert args[args.index("-contactMode") + 1] == "shadow-free"


def test_dropin_reads_cutoff_stack_scale_from_contact_scaling(tmp_path: Path) -> None:
    template = tmp_path / "temp.bifsif"
    template.mkdir()
    (template / "tmp.sif").write_text(
        """<?xml version='1.0'?>
<sif><settings>
<Groups><groupRatios contacts="1.2" dihedrals="1"/></Groups>
<Contacts method="cutoff" contactDistance="5.0">
<contactScaling name="stackingScale" scale="0.3"/>
</Contacts>
<dihedralNormalization dihedralCounting="1"/>
</settings></sif>
""",
        encoding="utf-8",
    )

    args = dropin.translate_smogcheck_args(["-i", "x.pdb", "-t", str(template)])

    assert args[args.index("-contactMode") + 1] == "cutoff"
    assert args[args.index("-contactStackScale") + 1] == "0.3"
    assert args[args.index("-dihedralCounting") + 1] == "1"


def test_dropin_translates_public_template_directory_files(tmp_path: Path) -> None:
    template = tmp_path / "SBM_AA+coulomb"
    template.mkdir()
    (template / "SBM_AA+coulomb.bif").write_text("<bif />\n", encoding="utf-8")
    (template / "SBM_AA+coulomb.nb").write_text("<nb />\n", encoding="utf-8")
    (template / "SBM_AA+coulomb.sif").write_text(
        """<?xml version='1.0'?>
<sif><settings>
<Contacts method="shadow" contactDistance="6" shadowRadius="1" shadowRadiusBonded="0.5"/>
</settings></sif>
""",
        encoding="utf-8",
    )

    args = dropin.translate_smogcheck_args(["-i", "x.pdb", "-t", str(template), "-opensmog"])

    assert "-OpenSMOG" in args
    assert "-AA" in args
    assert args[args.index("-templateBif") + 1] == str(template / "SBM_AA+coulomb.bif")
    assert args[args.index("-templateNb") + 1] == str(template / "SBM_AA+coulomb.nb")
    assert args[args.index("-contactMode") + 1] == "shadow"


def test_dropin_uses_tcg_template_as_final_model_template(tmp_path: Path) -> None:
    aa_template = tmp_path / "SBM_AA"
    ca_template = tmp_path / "SBM_CA+customNonbonded"
    aa_template.mkdir()
    ca_template.mkdir()
    (aa_template / "AA-whitford09.bif").write_text("<bif />\n", encoding="utf-8")
    (aa_template / "AA-whitford09.nb").write_text("<nb />\n", encoding="utf-8")
    (ca_template / "CA+customContacts.bif").write_text("<bif />\n", encoding="utf-8")
    (ca_template / "CA+customContacts.nb").write_text("<nb />\n", encoding="utf-8")

    args = dropin.translate_smogcheck_args(["-i", "x.pdb", "-t", str(aa_template), "-tCG", str(ca_template)])
    template_bif_indices = [idx for idx, arg in enumerate(args) if arg == "-templateBif"]

    assert "-CA" in args
    assert "-AA" not in args
    assert args[template_bif_indices[-1] + 1] == str(ca_template / "CA+customContacts.bif")


def test_opensmog_nonbond_settings_read_public_extras(tmp_path: Path) -> None:
    template = tmp_path / "Custom.bif"
    template.write_text("<bif />\n", encoding="utf-8")
    template.with_suffix(".sif").write_text(
        """<?xml version='1.0'?>
<sif><CustomNonBonded OpenSMOGparameters="epsilon,rmin,kappa"
OpenSMOGpotential="epsilon*exp(-(r-rmin)*kappa)" r_c="2.0"/></sif>
""",
        encoding="utf-8",
    )
    template.with_suffix(".extras").write_text(
        "nonbond_params < NB_1 NB_2 1 2.0 0.4 2.5\n",
        encoding="utf-8",
    )

    _constants, nonbond = smog2_native._opensmog_nonbond_settings(template)

    assert nonbond is not None
    assert nonbond["expression"] == "epsilon(type1,type2)*exp(-(r-rmin(type1,type2))*kappa(type1,type2))"
    assert nonbond["cutoff"] == "2.0"
    assert nonbond["parameter_rows"] == [("NB_1", "NB_2", ["2.0", "0.4", "2.5"])]
    assert smog2_native._opensmog_numeric_attr("2.5") == "2.50000e+00"


def test_dropin_detects_cutoff_gaussian_template(tmp_path: Path) -> None:
    template = tmp_path / "temp.bifsif"
    template.mkdir()
    (template / "tmp.sif").write_text(
        """<?xml version='1.0'?>
<sif>
<functions><function name="contact_gaussian" directive="pairs"/></functions>
<settings><Contacts method="cutoff" contactDistance="5.0"><contactScaling name="stackingScale" scale="0.3"/></Contacts></settings>
</sif>
""",
        encoding="utf-8",
    )

    args = dropin.translate_smogcheck_args(["-i", "x.pdb", "-t", str(template)])

    assert "-AAgaussian" in args
    assert args[args.index("-contactMode") + 1] == "cutoff-gaussian"


def test_dropin_translates_interactive_harness_choice() -> None:
    assert dropin.translate_smogcheck_args(["-i", "x.pdb"], stdin_text="0\n")[-1] == "-AA"
    assert dropin.translate_smogcheck_args(["-i", "x.pdb"], stdin_text="2\n")[-1] == "-CA"


def test_dropin_does_not_read_stdin_for_version_probe() -> None:
    assert dropin.translate_smogcheck_args(["-v"], stdin_text="share/settings/smog2.testlist\n") == ["-v"]


def test_dropin_main_forces_no_perl_fallback_and_writes_shadow_output(tmp_path, monkeypatch) -> None:
    contacts = tmp_path / "model.contacts"

    def fake_native_main(args: list[str]) -> int:
        assert "-SCMorig" not in args
        assert "-AA" in args
        contacts.write_text("1 2 1 5\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(dropin, "smog2_native_main", fake_native_main)
    monkeypatch.setenv("SMOG3_LEGACY_PERL_FALLBACK", "1")

    rc = dropin.main(["-SCMorig", "-i", "x.pdb", "-t", "SBM_AA", "-s", str(contacts)])

    assert rc == 0
    assert contacts.with_suffix(".contacts.ShadowOutput").read_text(encoding="utf-8") == "1 2 1 5\n"


def test_dropin_opensmog_v27_moves_custom_dihedrals_to_xml(tmp_path: Path) -> None:
    pdb = ROOT / "SMOG-CHECK" / "share" / "PDB.files" / "1F4N_v2.pdb"
    template = ROOT / "SMOG-CHECK" / "share" / "templates" / "SBM_AA+customContacts+customDihedrals" / "AA+customContacts+customDihedrals.bif"
    nb = template.with_suffix(".nb")
    atoms = smog2_native._parse_pdb_atoms(pdb)
    contacts = [(1, 47, tuple())]

    all_top = tmp_path / "all.top"
    dropin_top = tmp_path / "dropin.top"
    xml = tmp_path / "model.xml"
    smog2_native._write_case1_final_top(
        all_top,
        atoms,
        contacts,
        include_pairs=False,
        include_exclusions=False,
        template_path=template,
        nb_path=nb,
    )
    smog2_native._write_case1_final_top(
        dropin_top,
        atoms,
        contacts,
        include_pairs=False,
        include_exclusions=False,
        template_path=template,
        nb_path=nb,
        omit_proper_energy_groups={"sc_a"},
    )
    smog2_native._write_opensmog_xml(
        xml,
        atoms,
        contacts,
        "AA-CCD",
        template_path=template,
        nb_path=nb,
        generate_pair_exclusions=True,
    )

    all_func1 = sum(1 for row in _top_section_rows(all_top.read_text(encoding="utf-8"), "dihedrals") if row[4] == "1")
    dropin_func1 = sum(1 for row in _top_section_rows(dropin_top.read_text(encoding="utf-8"), "dihedrals") if row[4] == "1")
    custom_count = len(smog2_native._opensmog_custom_dihedral_items(atoms, template))
    assert all_func1 - dropin_func1 == 2 * custom_count
    assert "[ exclusions ]" not in dropin_top.read_text(encoding="utf-8")

    xml_text = xml.read_text(encoding="utf-8")
    assert '<exclusions generate="1"/>' in xml_text
    assert '<dihedrals_type name="dihedral_custom1">' in xml_text
    first_custom = re.search(r'<interaction i="[^"]+" j="[^"]+" k="[^"]+" l="[^"]+" theta0="([^"]+)" weight="([^"]+)" multiplicity="([^"]+)"/>', xml_text)
    assert first_custom is not None
    assert float(first_custom.group(2)) > 0.0


def test_dropin_dihe_defaults_and_type4_writer(tmp_path: Path) -> None:
    pdb = ROOT / "SMOG-CHECK" / "share" / "PDB.files" / "2ci2_v2.pdb"
    atoms = smog2_native._parse_pdb_atoms(pdb)
    top = tmp_path / "dihe4.top"

    smog2_native._write_case1_final_top(
        top,
        atoms,
        [(1, 20, tuple())],
        defaults_line="  1      1         no",
        proper_dihedral_func=4,
    )

    defaults = _top_section_rows(top.read_text(encoding="utf-8"), "defaults")
    dihedrals = _top_section_rows(top.read_text(encoding="utf-8"), "dihedrals")
    assert len(defaults[0]) == 3
    assert any(row[4] == "4" for row in dihedrals)
    assert all(row[4] != "1" for row in dihedrals)
