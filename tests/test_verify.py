from unittest.mock import Mock, patch

from smog3 import verify


def _cp(rc: int):
    m = Mock()
    m.returncode = rc
    return m


def test_verify_returns_nonzero_if_pytest_fails():
    with patch("smog3.verify.subprocess.run", side_effect=[_cp(1)]):
        assert verify.main([]) == 1


def test_verify_returns_zero_when_pytest_passes_and_full_not_required_without_perl_deps():
    with patch("smog3.verify.subprocess.run", side_effect=[_cp(0), _cp(1)]):
        assert verify.main([]) == 0


def test_verify_require_full_fails_without_perl_deps():
    with patch("smog3.verify.subprocess.run", side_effect=[_cp(0), _cp(1)]):
        assert verify.main(["--require-full"]) == 3
