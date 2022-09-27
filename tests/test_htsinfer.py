"""Unit tests for module ``htsinfer.py``"""

from os import linesep
from pathlib import Path

import pytest

from htsinfer.exceptions import (
    FileProblem,
    WorkEnvProblem,
)
from htsinfer.htsinfer import HtsInfer
from htsinfer.models import (
    ResultsLayout,
    ResultsOrientation,
    ResultsSource,
    ResultsStats,
    ResultsType,
    RunStates,
    StatesType,
    Results,
    Config,
    Args,
)
from tests.utils import (
    FILE_EMPTY,
    FILE_MATE_1,
    FILE_TRANSCRIPTS,
    RaiseFileProblem,
    RaiseMetadataWarning,
    RaiseOSError,
    CONFIG,
)


class TestHtsInfer:
    """Test ``HtsInfer`` class."""

    def test_init_no_args(self):
        """No arguments."""
        with pytest.raises(TypeError):
            HtsInfer()  # type: ignore

    def test_init_one_path(self):
        """One input file path."""
        CONFIG.args.path_1 = FILE_MATE_1
        test_instance = HtsInfer(config=CONFIG)
        assert test_instance.config.args.path_1 == FILE_MATE_1
        assert test_instance.state is RunStates.OKAY

    def test_init_two_paths(self):
        """Two input file paths."""
        CONFIG.args.path_2 = FILE_MATE_1
        test_instance = HtsInfer(config=CONFIG)
        assert test_instance.config.args.path_1 == FILE_MATE_1
        assert test_instance.config.args.path_2 == FILE_MATE_1
        assert test_instance.state is RunStates.OKAY

    def test_evaluate(self, tmpdir, monkeypatch):
        """No warnings."""
        CONFIG.args.path_2 = None
        CONFIG.args.tmp_dir = Path(tmpdir)
        CONFIG.args.out_dir = Path(tmpdir)
        test_instance = HtsInfer(config=CONFIG)
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.evaluate',
            ResultsSource,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_stats.GetLibStats.evaluate',
            ResultsStats,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_type.GetLibType.evaluate',
            ResultsType,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_layout.GetReadLayout.evaluate',
            ResultsLayout,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_orientation.GetOrientation.evaluate',
            ResultsOrientation,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.OKAY

    def test_evaluate_lib_type_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in library type determination."""
        test_instance = HtsInfer(config=CONFIG)
        monkeypatch.setattr(
            'htsinfer.get_library_type.GetLibType.evaluate',
            RaiseMetadataWarning,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.evaluate',
            ResultsSource,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_stats.GetLibStats.evaluate',
            ResultsStats,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_layout.GetReadLayout.evaluate',
            ResultsLayout,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_orientation.GetOrientation.evaluate',
            ResultsOrientation,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_read_orient_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in read orientation determination."""
        test_instance = HtsInfer(config=CONFIG)
        monkeypatch.setattr(
            'htsinfer.get_read_orientation.GetOrientation.evaluate',
            RaiseMetadataWarning,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.evaluate',
            ResultsSource,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_stats.GetLibStats.evaluate',
            ResultsStats,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_type.GetLibType.evaluate',
            ResultsType,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_layout.GetReadLayout.evaluate',
            ResultsLayout,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_read_layout_metadata_warning(self, monkeypatch, tmpdir):
        """Metadata warning in read layout determination."""
        test_instance = HtsInfer(config=CONFIG)
        monkeypatch.setattr(
            'htsinfer.get_read_layout.GetReadLayout.evaluate',
            RaiseMetadataWarning,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.evaluate',
            ResultsSource,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_stats.GetLibStats.evaluate',
            ResultsStats,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_type.GetLibType.evaluate',
            ResultsType,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_orientation.GetOrientation.evaluate',
            ResultsOrientation,
        )
        test_instance.evaluate()
        assert test_instance.state is RunStates.WARNING

    def test_evaluate_file_problem(self, tmpdir):
        """File problem due to empty file."""
        CONFIG.args.path_1 = FILE_EMPTY
        test_instance = HtsInfer(config=CONFIG)
        test_instance.evaluate()
        assert test_instance.state is RunStates.ERROR

    def test_evaluate_work_env_problem(self, tmpdir):
        """Cannot create work environment."""
        test_instance = HtsInfer(config=CONFIG)
        test_instance.out_dir = Path(".")
        test_instance.evaluate()
        assert test_instance.state is RunStates.ERROR

    def test_prepare_env_default(self, tmpdir):
        """Test default behavior."""
        CONFIG.args.path_1 = FILE_MATE_1
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        assert Path(test_instance.config.args.out_dir).is_dir()
        assert Path(test_instance.config.args.tmp_dir).is_dir()

    def test_prepare_env_out_dir_creation_fails(self):
        """Creation of output directory fails."""
        test_instance = HtsInfer(config=CONFIG)
        test_instance.config.args.out_dir = Path(".")
        with pytest.raises(WorkEnvProblem):
            test_instance.prepare_env()

    def test_prepare_env_tmp_dir_creation_fails(self, tmpdir):
        """Creation of temporary directory fails."""
        CONFIG.args.out_dir = tmpdir
        test_instance = HtsInfer(config=CONFIG)
        test_instance.config.args.tmp_dir = Path(".")
        with pytest.raises(WorkEnvProblem):
            test_instance.prepare_env()

    def test_process_inputs_default(self, tmpdir):
        """Test default behavior."""
        CONFIG.args.path_1 = FILE_MATE_1
        CONFIG.args.path_2 = FILE_MATE_1
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.tmp_dir = tmpdir
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        test_instance.process_inputs()
        test_instance.clean_up()

    def test_process_input_param_orient_min_fraction_value_error(self, tmpdir):
        """Invalid value for read orientation min fraction parameter."""
        CONFIG.args.read_orientation_min_fraction = 0.49
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        with pytest.raises(ValueError):
            test_instance.process_inputs()
        test_instance.clean_up()

    def test_process_inputs_file_problem_empty(self, tmpdir):
        """File validation fails because input file is empty."""
        CONFIG.args.path_1 = FILE_EMPTY
        CONFIG.args.path_2 = None
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.read_orientation_min_fraction = 0.75
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        with pytest.raises(FileProblem):
            test_instance.process_inputs()
        test_instance.clean_up()

    def test_process_transcripts_file_unzipped(self, tmpdir):
        """Transcripts file is not compressed."""
        CONFIG.args.path_1 = FILE_MATE_1
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.transcripts_file = FILE_TRANSCRIPTS
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        test_instance.process_inputs()
        test_instance.clean_up()

    def test_process_transcripts_file_problem_empty(self, monkeypatch, tmpdir):
        """File validation fails because transcripts file is empty."""
        CONFIG.args.transcripts_file = FILE_EMPTY
        test_instance = HtsInfer(config=CONFIG)
        monkeypatch.setattr(
            'shutil.copyfileobj',
            RaiseFileProblem,
        )
        test_instance.prepare_env()
        with pytest.raises(FileProblem):
            test_instance.process_inputs()
        test_instance.clean_up()

    def test_get_library_stats_default(self):
        """Test default behavior."""
        CONFIG.args.path_1 = FILE_MATE_1
        test_instance = HtsInfer(config=CONFIG)
        test_instance.get_library_stats()
        assert (
            test_instance.config.results.library_stats.file_1.read_length.min == 150
        )
        assert (
            test_instance.config.results.library_stats.file_1.read_length.max == 150
        )

    def test_get_library_type_default(self):
        """Test default behavior."""
        test_instance = HtsInfer(config=CONFIG)
        test_instance.get_library_type()
        assert (
            test_instance.config.results.library_type.file_1 is
            StatesType.first_mate
        )
        assert test_instance.config.results.library_type is not None

    def test_clean_up_keep_none(self, tmpdir):
        """Remove all data."""
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        test_instance.state = RunStates.OKAY
        test_instance.clean_up()
        assert not Path(test_instance.config.args.out_dir).is_dir()
        assert not Path(test_instance.config.args.tmp_dir).is_dir()

    def test_clean_up_keep_results(self, tmpdir):
        """Remove temporary data."""
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        test_instance.state = RunStates.WARNING
        test_instance.clean_up()
        assert Path(test_instance.config.args.out_dir).is_dir()
        assert not Path(test_instance.config.args.tmp_dir).is_dir()

    def test_clean_up_keep_all(self, tmpdir):
        """Remove no data."""
        test_instance = HtsInfer(config=CONFIG)
        test_instance.prepare_env()
        test_instance.state = RunStates.ERROR
        test_instance.clean_up()
        assert Path(test_instance.config.args.out_dir).is_dir()
        assert Path(test_instance.config.args.tmp_dir).is_dir()

    def test_clean_up_out_dir_removal_fails(self, monkeypatch, tmpdir):
        """Output directory removal fails."""
        test_instance = HtsInfer(config=CONFIG)
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
        test_instance = HtsInfer(config=CONFIG)
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
        arguments = Args(path_1=FILE_MATE_1)
        results = Results()
        configs = Config(
            args=arguments,
            results=results,
        )
        test_instance = HtsInfer(config=configs)
        test_instance.print()
        captured = capsys.readouterr()
        assert captured.out == ('''{
   "library_source": {
      "file_1": {
         "short_name": null,
         "taxon_id": null
      },
      "file_2": {
         "short_name": null,
         "taxon_id": null
      }
   },
   "library_stats": {
      "file_1": {
         "read_length": {
            "max": null,
            "min": null
         }
      },
      "file_2": {
         "read_length": {
            "max": null,
            "min": null
         }
      }
   },
   "library_type": {
      "file_1": null,
      "file_2": null,
      "relationship": null
   },
   "read_layout": {
      "file_1": {
         "adapt_3": null
      },
      "file_2": {
         "adapt_3": null
      }
   },
   "read_orientation": {
      "file_1": null,
      "file_2": null,
      "relationship": null
   }
}''') + linesep
        assert captured.err == ""
