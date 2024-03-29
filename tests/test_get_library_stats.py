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
    RaiseError,
    CONFIG,
)


class TestGetOrientation:
    """Test ``GetOrientation`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        test_instance = GetLibStats(config=CONFIG)
        assert test_instance.paths[0] == FILE_MATE_1

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        CONFIG.args.path_2_processed = FILE_MATE_2
        test_instance = GetLibStats(config=CONFIG)
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2

    def test_evaluate_single(self, tmpdir):
        """Get library statistics for a single-end library."""
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.path_2_processed = None
        test_instance = GetLibStats(
            config=CONFIG,
        )
        results = test_instance.evaluate()
        assert results == ResultsStats(
            file_1=Stats(read_length=ReadLength(
                min=150, max=150, mean=150.0, median=150, mode=150
            )),
            file_2=Stats(read_length=ReadLength(
                min=None, max=None, mean=None, median=None, mode=None
            )),
        )

    def test_evaluate_paired(self, tmpdir):
        """Get library statistics for a paired-end library."""
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.path_2_processed = FILE_MATE_2
        test_instance = GetLibStats(
            config=CONFIG,
        )
        results = test_instance.evaluate()
        assert results == ResultsStats(
            file_1=Stats(read_length=ReadLength(
                min=150, max=150, mean=150.0, median=150, mode=150
            )),
            file_2=Stats(read_length=ReadLength(
                min=150, max=150, mean=150.0, median=150, mode=150
            )),
        )

    def test_fastq_get_stats_read_length(self):
        """Get min/max read length of a FASTQ file."""
        results = GetLibStats.fastq_get_stats_read_length(fastq=FILE_MATE_1)
        assert results == (150, 150, 150.0, 150, 150)

    def test_fastq_get_stats_read_length_file_problem(self, monkeypatch):
        """File problem when accessing FASTQ file."""
        monkeypatch.setattr(
            'Bio.SeqIO.parse',
            RaiseError(exc=OSError),
        )
        with pytest.raises(FileProblem):
            GetLibStats.fastq_get_stats_read_length(fastq=FILE_MATE_1)
