"""Unit tests for module ``htsinfer.py``"""

from pathlib import Path

import pytest

from htsinfer.exceptions import (
    FileProblem,
    WorkEnvProblem,
)
from htsinfer.htsinfer import HtsInfer
from htsinfer.models import (
    RunStates,
    StatesType,
)
from tests.utils import (
    FILE_EMPTY,
    FILE_MATE_1,
    RaiseMetadataWarning,
    RaiseOSError,
)


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
        assert test_instance.state is RunStates.OKAY

    def test_init_two_paths(self):
        """Two input file paths."""
        test_instance = HtsInfer(path_1=FILE_MATE_1, path_2=FILE_MATE_1)
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.path_2 == FILE_MATE_1
        assert test_instance.state is RunStates.OKAY

    def test_evaluate(self, tmpdir):
        """No warnings."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.results.library_type is not None
        assert test_instance.state is RunStates.OKAY

    def test_evaluate_lib_type_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in library type determination."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_library_type',
            RaiseMetadataWarning,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_lib_source_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in library source determination."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_library_source',
            RaiseMetadataWarning,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_read_orient_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in read orientation determination."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_read_orientation',
            RaiseMetadataWarning,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_read_layout_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in read layout determination."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_read_layout',
            RaiseMetadataWarning,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_file_problem(self, tmpdir):
        """File problem due to empty file."""
        test_instance = HtsInfer(
            path_1=FILE_EMPTY,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.ERROR

    def test_evaluate_work_env_problem(self, tmpdir):
        """Cannot create work environment."""
        test_instance = HtsInfer(path_1=FILE_EMPTY)
        test_instance.out_dir = Path(".")
        test_instance.evaluate()
        assert test_instance.state is RunStates.ERROR

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

    def test_process_inputs_default(self, tmpdir):
        """Test default behavior."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            path_2=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.process_inputs()
        test_instance.clean_up()

    def test_process_inputs_file_problem_empty(self, tmpdir):
        """File validation fails because input file is empty."""
        test_instance = HtsInfer(
            path_1=FILE_EMPTY,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        with pytest.raises(FileProblem):
            test_instance.process_inputs()
        test_instance.clean_up()

    def test_get_library_type_default(self):
        """Test default behavior."""
        test_instance = HtsInfer(path_1=FILE_MATE_1)
        test_instance.get_library_type()
        assert (
            test_instance.results.library_type.file_1 is
            StatesType.first_mate
        )
        assert test_instance.results.library_type is not None

    def test_clean_up_keep_none(self, tmpdir):
        """Remove all data."""
        test_instance = HtsInfer(
            path_1=FILE_MATE_1,
            out_dir=tmpdir,
            tmp_dir=tmpdir,
        )
        test_instance.prepare_env()
        test_instance.state = RunStates.OKAY
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
        test_instance.state = RunStates.WARNING
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
        test_instance.state = RunStates.ERROR
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
        test_instance.state = RunStates.OKAY
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
        test_instance.state = RunStates.WARNING
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
            '"library_type": {'
            '"file_1": null, "file_2": null, "relationship": null'
            '}, '
            '"library_source": {}, '
            '"read_orientation": {}, '
            '"read_layout": {'
            '"file_1": '
            '{'
            '"adapt_3": null'
            '}, '
            '"file_2": '
            '{'
            '"adapt_3": null'
            '}'
            '}'
            '}'
        )
        assert captured.err == ""
