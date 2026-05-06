import subprocess
from pathlib import Path


def test_install_compat_generates_smogcheck_friendly_wrappers(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[1]
    bin_dir = tmp_path / "bin"

    subprocess.run(
        [
            "python",
            "-c",
            "from smog3.install_compat import main; raise SystemExit(main())",
            "--bin-dir",
            str(bin_dir),
            "--repo-root",
            str(repo_root),
            "--perl",
            "perl",
        ],
        check=True,
        env={**__import__("os").environ, "PYTHONPATH": "src"},
        cwd=repo_root,
    )

    smog2 = (bin_dir / "smog2").read_text()
    assert "export PERLLIB=" in smog2
    assert "export SMOG_PATH=" in smog2
    assert "src/smogv2 \"$@\"" in smog2

    tail = smog2.strip().splitlines()[-1].split()
    assert tail[0] == "perl"
    assert tail[1].endswith("src/smogv2")

    assert (bin_dir / "smog_tablegen").exists()
    assert (bin_dir / "smog_scale-energies").exists()
