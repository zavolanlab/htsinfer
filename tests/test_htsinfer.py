"""Unit tests for module ``htsinfer.py``"""

from pathlib import Path

import pytest

from htsinfer.exceptions import WorkEnvProblem
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
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.results.library_type is None
        assert test_instance.state is RunStates.okay

    def test_evaluate_work_env_problem(self, tmpdir):
        """Cannot create work environment."""
        test_instance = HtsInfer(path_1=FILE_EMPTY)
        test_instance.out_dir = Path(".")
        test_instance.evaluate()
        assert test_instance.state is RunStates.error

    def test_prepare_env_default(self, tmpdir):
        """Test default behavior."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        assert Path(test_instance.out_dir).is_dir()
        assert Path(test_instance.tmp_dir).is_dir()

    def test_prepare_env_out_dir_creation_fails(self):
        """Creation of output directory fails."""
        test_instance = HtsInfer(path_1=FILE_MATE_1)
        test_instance.out_dir = Path(".")
        with pytest.raises(WorkEnvProblem):
            test_instance.prepare_env()

    def test_prepare_env_tmp_dir_creation_fails(self, tmpdir):
        """Creation of temporary directory fails."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
        )
        test_instance.tmp_dir = Path(".")
        with pytest.raises(WorkEnvProblem):
            test_instance.prepare_env()

    def test_clean_up_keep_none(self, tmpdir):
        """Remove all data."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.state = RunStates.okay
        test_instance.clean_up()
        assert not Path(test_instance.out_dir).is_dir()
        assert not Path(test_instance.tmp_dir).is_dir()

    def test_clean_up_keep_results(self, tmpdir):
        """Remove temporary data."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.state = RunStates.warning
        test_instance.clean_up()
        assert Path(test_instance.out_dir).is_dir()
        assert not Path(test_instance.tmp_dir).is_dir()

    def test_clean_up_keep_all(self, tmpdir):
        """Remove no data."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.state = RunStates.error
        test_instance.clean_up()
        assert Path(test_instance.out_dir).is_dir()
        assert Path(test_instance.tmp_dir).is_dir()

    def test_clean_up_out_dir_removal_fails(self, monkeypatch, tmpdir):
        """Output directory removal fails."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.state = RunStates.okay
        monkeypatch.setattr(
            'shutil.rmtree',
            RaiseOSError,
        )
        with pytest.raises(WorkEnvProblem):
            test_instance.clean_up()

    def test_clean_up_tmp_dir_removal_fails(self, monkeypatch, tmpdir):
        """Temporary directory removal fails."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.state = RunStates.warning
        monkeypatch.setattr(
            'shutil.rmtree',
            RaiseOSError,
        )
        with pytest.raises(WorkEnvProblem):
            test_instance.clean_up()

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
