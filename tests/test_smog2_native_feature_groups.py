from pathlib import Path
import xml.etree.ElementTree as ET

import smog3.cli as cli


OPENSMOG_CASES = [
    ("1AKEapo_v2", ["-AA"]),
    ("DNA.terminal.BMG", ["-AA"]),
    ("1A01-AMP", ["-AAgaussian"]),
    ("2ci2_v2", ["-CA"]),
    ("2ci2_v2", ["-CAgaussian"]),
    ("1AKEapo_v2", ["-AA"]),  # shadow-free variant approximated by AA OpenSMOG path in native slice
    ("2ci2_v2", ["-AACC1"]),
    ("tRNA", ["-AACC1"]),
    ("terminaltest.BOND", ["-AACC1"]),
    ("1F4N_v2", ["-AACC1"]),
    ("1F4N_v2", ["-AACCD"]),
    ("2ci2_v2.freecoor", ["-AA"]),
    ("RNA+protein", ["-AADIHE"]),
    ("RNA+protein", ["-AADIHE4"]),
]


def _count_atoms(pdb: Path) -> int:
    return sum(1 for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM")))


def test_g96_and_opensmog_groups_run_natively_without_perl(monkeypatch, tmp_path: Path):
    pdb_root = Path(__file__).resolve().parents[1] / "SMOG-CHECK" / "share" / "PDB.files"

    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked: {argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)

    # Group 1: G96 output writer
    g96_pdb = pdb_root / "2ci2_v2.pdb"
    g96_top = tmp_path / "g96_case.top"
    g96_coord = tmp_path / "g96_case.g96"
    g96_ndx = tmp_path / "g96_case.ndx"
    g96_contacts = tmp_path / "g96_case.contacts"
    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["smog2", "-i", str(g96_pdb), "-AA", "-g96", "-o", str(g96_top), "-g", str(g96_coord), "-n", str(g96_ndx), "-s", str(g96_contacts)],
    )
    assert cli.smog2_main() == 0
    g96_lines = g96_coord.read_text().splitlines()
    natoms = _count_atoms(g96_pdb)
    assert g96_lines[0] == "TITLE"
    assert "POSITION" in g96_lines
    assert "BOX" in g96_lines
    pos_idx = g96_lines.index("POSITION")
    end_idx = g96_lines.index("END", pos_idx + 1)
    assert (end_idx - pos_idx - 1) == natoms

    # Group 2: OpenSMOG XML generation (cases 94-104, 106, 111, 113)
    for idx, (stem, mode_args) in enumerate(OPENSMOG_CASES, start=1):
        pdb = pdb_root / f"{stem}.pdb"
        base = tmp_path / f"opensmog_{idx}"
        top = base.with_suffix(".top")
        gro = base.with_suffix(".gro")
        ndx = base.with_suffix(".ndx")
        contacts = base.with_suffix(".contacts")
        xml = base.with_suffix(".xml")
        monkeypatch.setattr(
            cli.sys,
            "argv",
            ["smog2", "-i", str(pdb), *mode_args, "-OpenSMOG", "-OpenSMOGxml", str(xml), "-o", str(top), "-g", str(gro), "-n", str(ndx), "-s", str(contacts)],
        )
        assert cli.smog2_main() == 0

        natoms = _count_atoms(pdb)
        top_text = top.read_text()
        assert "[ atoms ]" in top_text and "[ bonds ]" in top_text and "[ molecules ]" in top_text
        assert int(gro.read_text().splitlines()[1].strip()) == natoms
        assert len([x for x in ndx.read_text().split() if x.isdigit()]) == natoms

        tree = ET.parse(xml)
        root = tree.getroot()
        assert root.tag == "OpenSMOG"
        ff = root.find("ForceField")
        assert ff is not None
        particles = ff.find("Particles")
        assert particles is not None
        particle_nodes = particles.findall("Particle")
        assert int(particles.attrib["count"]) == natoms == len(particle_nodes)
        contacts_node = ff.find("Contacts")
        assert contacts_node is not None
        contact_nodes = contacts_node.findall("Contact")
        file_contacts = [ln for ln in contacts.read_text().splitlines() if ln.strip()]
        assert int(contacts_node.attrib["count"]) == len(contact_nodes) == len(file_contacts)

    assert called["perl"] is False
