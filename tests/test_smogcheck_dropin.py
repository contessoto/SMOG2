from pathlib import Path

from smog3 import smogcheck_dropin_smog2 as dropin


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
