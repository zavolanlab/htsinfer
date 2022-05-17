"""Unit tests for module ``get_read_library_stats.py``."""

import pytest

from htsinfer.exceptions import FileProblem
from htsinfer.get_library_stats import GetLibStats
from htsinfer.models import (
    ReadLength,
    ResultsStats,
    Stats,
)
from tests.utils import (
    FILE_MATE_1,
    FILE_MATE_2,
    RaiseOSError,
)


class TestGetOrientation:
    """Test ``GetOrientation`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        test_instance = GetLibStats(paths=(FILE_MATE_1, None))
        assert test_instance.paths[0] == FILE_MATE_1

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        test_instance = GetLibStats(paths=(FILE_MATE_1, FILE_MATE_2))
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2

    def test_evaluate_single(self, tmpdir):
        """Get library statistics for a single-end library."""
        test_instance = GetLibStats(
            paths=(FILE_MATE_1, None),
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsStats(
            file_1=Stats(read_length=ReadLength(min=150, max=150)),
            file_2=Stats(read_length=ReadLength(min=None, max=None)),
        )

    def test_evaluate_paired(self, tmpdir):
        """Get library statistics for a paired-end library."""
        test_instance = GetLibStats(
            paths=(FILE_MATE_1, FILE_MATE_2),
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsStats(
            file_1=Stats(read_length=ReadLength(min=150, max=150)),
            file_2=Stats(read_length=ReadLength(min=150, max=150)),
        )

    def test_fastq_get_min_max_read_length(self):
        """Get min/max read length of a FASTQ file."""
        results = GetLibStats.fastq_get_min_max_read_length(fastq=FILE_MATE_1)
        assert results == (150, 150)

    def test_fastq_get_min_max_read_length_file_problem(self, monkeypatch):
        """File problem when accessing FASTQ file."""
        monkeypatch.setattr(
            'Bio.SeqIO.parse',
            RaiseOSError,
        )
        with pytest.raises(FileProblem):
            GetLibStats.fastq_get_min_max_read_length(fastq=FILE_MATE_1)
