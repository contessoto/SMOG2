from pathlib import Path

import smog3.cli as cli


def _count_section(top_text: str, section: str) -> int:
    in_section = False
    n = 0
    for ln in top_text.splitlines():
        s = ln.strip()
        if s.startswith("[") and s.endswith("]"):
            in_section = s == f"[ {section} ]"
            continue
        if in_section and s and not s.startswith(";"):
            n += 1
    return n


def test_gaussian_batch_cases_use_native_with_explicit_output_options(monkeypatch, tmp_path: Path):
    # smog2.testlist cases 21-25: AA default-gaussian
    stems = ["1A01-AMP", "1A01-ADP", "1A01-ATP", "1A01-ADP+RNA", "2FP4-GDP"]
    pdb_root = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files"

    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked: {argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)

    for stem in stems:
        pdb = pdb_root / f"{stem}.pdb"
        out_base = tmp_path / stem.replace("+", "_")
        top = out_base.with_suffix(".top")
        gro = out_base.with_suffix(".gro")
        ndx = out_base.with_suffix(".ndx")
        contacts = out_base.with_suffix(".contacts")

        monkeypatch.setattr(
            cli.sys,
            "argv",
            ["smog2", "-i", str(pdb), "-AAgaussian", "-o", str(top), "-g", str(gro), "-n", str(ndx), "-s", str(contacts)],
        )
        rc = cli.smog2_main()
        assert rc == 0

        top_text = top.read_text()
        gro_lines = gro.read_text().splitlines()
        ndx_ids = [x for x in ndx.read_text().split() if x.isdigit()]
        contact_lines = [ln for ln in contacts.read_text().splitlines() if ln.strip()]
        pdb_atoms = sum(1 for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM")))

        assert int(gro_lines[1].strip()) == pdb_atoms
        assert _count_section(top_text, "atoms") == pdb_atoms
        assert _count_section(top_text, "bonds") == max(0, pdb_atoms - 1)
        assert _count_section(top_text, "molecules") == 1
        assert len(ndx_ids) == pdb_atoms
        assert len(contact_lines) == max(0, pdb_atoms - 3)

        assert contacts.name.endswith(".contacts")
        assert gro.name.endswith(".gro")
        assert top.name.endswith(".top")
        assert ndx.name.endswith(".ndx")

    assert called["perl"] is False
