from pathlib import Path

from smog3.parity_report import compare_digest_sets, digest_files


def test_digest_and_compare(tmp_path: Path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir(); b.mkdir()

    (a / "one.txt").write_text("abc")
    (a / "same.txt").write_text("same")
    (b / "same.txt").write_text("same")
    (b / "one.txt").write_text("xyz")
    (b / "only_b.txt").write_text("b")

    da = digest_files(a, ["*.txt"])
    db = digest_files(b, ["*.txt"])
    rep = compare_digest_sets(da, db)

    assert rep["match"] is False
    assert rep["only_in_a"] == []
    assert rep["only_in_b"] == ["only_b.txt"]
    assert rep["mismatched"][0]["path"] == "one.txt"
