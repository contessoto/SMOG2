from pathlib import Path

import smog3.cli as cli


MODES = {"shadow", "shadow-free", "cutoff", "cutoff-gaussian"}


def _parse_atoms(pdb: Path):
    atoms = []
    for ln in pdb.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            atoms.append((x, y, z))
    return atoms


def _expected_contacts(atoms, cutoff: float, min_sep: int):
    out = []
    c2 = cutoff * cutoff
    n = len(atoms)
    step = max(1, n // 1200)
    indices = list(range(0, n, step))
    for ii, i in enumerate(indices):
        xi, yi, zi = atoms[i]
        for j in indices[ii + 1:]:
            if (j - i) <= min_sep:
                continue
            dx = xi - atoms[j][0]
            dy = yi - atoms[j][1]
            dz = zi - atoms[j][2]
            if dx * dx + dy * dy + dz * dz <= c2:
                out.append((str(i + 1), str(j + 1)))
    return out


def _iter_cases_56_92(testlist: Path):
    for raw in testlist.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if ";" not in line:
            continue
        left, cid_txt = line.split(";", 1)
        try:
            cid = int(cid_txt.strip())
        except ValueError:
            continue
        if cid < 56 or cid > 92:
            continue
        cols = left.split()
        stem, model, mode = cols[0], cols[1], cols[2]
        if mode not in MODES:
            continue
        param = float(cols[3])
        yield cid, stem, model, mode, param


def test_shadow_cutoff_case_group_runs_native_without_perl(monkeypatch, tmp_path: Path):
    called = {"perl": False}

    def fail_perl(argv):
        called["perl"] = True
        raise AssertionError(f"Perl fallback invoked: {argv}")

    monkeypatch.setattr(cli, "run_smog2", fail_perl)

    root = Path(__file__).resolve().parents[1]
    testlist = root / "SMOG-CHECK" / "share" / "settings" / "smog2.testlist"
    pdb_root = root / "SMOG-CHECK" / "share" / "PDB.files"

    covered = []
    for cid, stem, model, mode, param in _iter_cases_56_92(testlist):
        flag = "-AA" if model == "AA" else "-CA"
        pdb = pdb_root / f"{stem}.pdb"
        base = tmp_path / f"c{cid}_{stem.replace('.', '_')}"
        top = base.with_suffix('.top')
        gro = base.with_suffix('.gro')
        ndx = base.with_suffix('.ndx')
        contacts = base.with_suffix('.contacts')

        monkeypatch.setattr(
            cli.sys,
            "argv",
            ["smog2", "-i", str(pdb), flag, "-contactMode", mode, "-contactParam", str(param), "-o", str(top), "-g", str(gro), "-n", str(ndx), "-s", str(contacts)],
        )
        rc = cli.smog2_main()
        assert rc == 0

        natoms = sum(1 for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM", "HETATM")))
        assert int(gro.read_text().splitlines()[1].strip()) == natoms
        assert len([x for x in ndx.read_text().split() if x.isdigit()]) == natoms
        t = top.read_text()
        assert "[ atoms ]" in t and "[ bonds ]" in t and "[ molecules ]" in t

        rows = [tuple(ln.split()[:2]) for ln in contacts.read_text().splitlines() if ln.strip()]
        atoms_xyz = _parse_atoms(pdb)
        min_sep = 3 if model == "AA" else 2
        expected = _expected_contacts(atoms_xyz, abs(param), min_sep)
        assert rows == expected
        covered.append(cid)

    assert covered == list(range(56, 93))
    assert called["perl"] is False
