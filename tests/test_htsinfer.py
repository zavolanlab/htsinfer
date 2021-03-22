"""Unit tests for module ``htsinfer.py``"""

from pathlib import Path

import pytest

from htsinfer.htsinfer import (HtsInfer, RunStates)

# Test parameters
TEST_FILES_DIR = Path(__file__).resolve().parent / "files"
FILE_EMPTY = TEST_FILES_DIR / "empty.fastq"
FILE_MATE_1 = TEST_FILES_DIR / "first_mate.fastq"


class RaiseOSError:
    """Raise ``OSError``."""
    def __init__(self, *args, **kwargs):
        raise OSError


class TestHtsInfer:
    """Test ``HtsInfer`` class."""

    def test_init_no_args(self):
        """No arguments."""
        with pytest.raises(TypeError):
            HtsInfer()  # type: ignore

    def test_init_one_path(self):
        """One input file path."""
        test_instance = HtsInfer(path_1=FILE_MATE_1)
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.state is RunStates.okay

    def test_init_two_paths(self):
        """Two input file paths."""
        test_instance = HtsInfer(path_1=FILE_MATE_1, path_2=FILE_MATE_1)
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.path_2 == FILE_MATE_1
        assert test_instance.state is RunStates.okay

    def test_evaluate(self, tmpdir):
        """No warnings."""
        test_instance = HtsInfer(path_1=FILE_MATE_1)
        test_instance.evaluate()
        assert test_instance.results.library_type is None
        assert test_instance.state is RunStates.okay

    def test_print(self, capsys):
        """Test default behavior."""
        test_instance = HtsInfer(path_1=FILE_MATE_1)
        test_instance.print()
        captured = capsys.readouterr()
        assert captured.out == (
            '{'
            '"library_type": null, '
            '"library_source": null, '
            '"read_orientation": null, '
            '"read_layout": null'
            '}'
        )
        assert captured.err == ""
