import sys
import pytest
from seismic_forward.cli import main, get_version


def test_get_version():
    """Test that get_version returns a version string."""
    version = get_version()
    assert isinstance(version, str)
    assert len(version) > 0


def test_version_flag(monkeypatch):
    """Test that --version flag outputs version information and exits."""
    monkeypatch.setattr(sys, "argv", ["seismic_forward", "--version"])

    with pytest.raises(SystemExit) as exc_info:
        main()

    # argparse's version action exits with code 0
    assert exc_info.value.code == 0


def test_help_flag(monkeypatch, capsys):
    """Test that --help flag outputs help information and exits."""
    monkeypatch.setattr(sys, "argv", ["seismic_forward", "--help"])

    with pytest.raises(SystemExit) as exc_info:
        main()

    # argparse's help action exits with code 0
    assert exc_info.value.code == 0

    # Check that help output contains expected information
    captured = capsys.readouterr()
    output = (captured.out + captured.err).lower()
    assert "seismic_forward" in output or "usage" in output


def test_no_arguments(monkeypatch, capsys):
    """Test that running without arguments shows help and error message."""
    monkeypatch.setattr(sys, "argv", ["seismic_forward"])

    with pytest.raises(SystemExit) as exc_info:
        main()

    # Should exit with error code 1
    assert exc_info.value.code == 1

    # Check that help and error message is shown
    captured = capsys.readouterr()
    output = captured.out + captured.err
    assert "modelfile must be provided" in output.lower()


def test_nonexistent_modelfile(monkeypatch):
    """Test that providing a nonexistent file results in an error."""
    monkeypatch.setattr(sys, "argv", ["seismic_forward", "nonexistent_file.xml"])

    with pytest.raises(SystemExit) as exc_info:
        main()

    # Should exit with error code 1
    assert exc_info.value.code == 1
