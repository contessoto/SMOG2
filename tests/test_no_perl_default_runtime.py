import smog3.cli as cli


def test_smog2_unsupported_fails_loudly_without_perl(monkeypatch, capsys):
    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError("Perl fallback should not be called by default")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)
    monkeypatch.setattr(cli.sys, "argv", ["smog2", "-i", "missing.pdb", "-AA", "-unknownFlag"])
    rc = cli.smog2_main()
    out = capsys.readouterr().out
    assert rc == 2
    assert "Unsupported smog2 options for native runtime" in out
    assert called["perl"] is False


def test_helper_commands_do_not_use_perl_by_default(monkeypatch, tmp_path):
    called = {"perl": False, "tool": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError("Perl fallback should not be called")

    def fail_tool(tool_name, argv):
        called["tool"] = True
        raise AssertionError("Perl tool wrapper should not be called")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)
    monkeypatch.setattr(cli, "run_tool", fail_tool)

    for fn, name in [
        (cli.smog_tablegen_main, "smog_tablegen"),
        (cli.smog_scale_energies_main, "smog_scale-energies"),
        (cli.smog_extract_main, "smog_extract"),
        (cli.smog_adjustpdb_main, "smog_adjustPDB"),
        (cli.smog_editgro_main, "smog_editgro"),
        (cli.smog_ions_main, "smog_ions"),
        (cli.smog_modifyxml_main, "smog_modifyXML"),
    ]:
        monkeypatch.setattr(cli.sys, "argv", [name])
        try:
            fn()
        except (SystemExit, FileNotFoundError):
            pass
    assert called["perl"] is False

