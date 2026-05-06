from smog3 import compat


def test_root_layout_expected():
    assert (compat.ROOT / "src" / "smogv2").exists()
    assert (compat.ROOT / "src" / "tools" / "tablegen").exists()
