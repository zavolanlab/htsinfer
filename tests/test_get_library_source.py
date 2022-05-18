"""Unit tests for module ``get_read_library_source.py``."""

from htsinfer.get_library_source import GetLibSource
from htsinfer.models import (
    ResultsSource,
    Source,
)
from tests.utils import (
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_TRANSCRIPTS,
    SOURCE_HUMAN,
)


class TestGetOrientation:
    """Test ``GetOrientation`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
        )
        assert test_instance.paths[0] == FILE_MATE_1

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, FILE_MATE_2),
            transcripts_file=FILE_TRANSCRIPTS,
        )
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2

    def test_evaluate_single(self, monkeypatch, tmpdir):
        """Get library statistics for a single-end library."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_source',
            lambda *args, **kwargs: SOURCE_HUMAN,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=SOURCE_HUMAN,
            file_2=Source(),
        )

    def test_evaluate_paired(self, monkeypatch, tmpdir):
        """Get library statistics for a paired-end library."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, FILE_MATE_2),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_source',
            lambda *args, **kwargs: SOURCE_HUMAN,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=SOURCE_HUMAN,
            file_2=SOURCE_HUMAN,
        )
