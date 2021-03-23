"""Unit tests for module ``subset_fastq.py``."""

import pytest

from htsinfer.exceptions import FileProblem
from htsinfer.subset_fastq import SubsetFastq
from tests.utils import (
    FILE_EMPTY,
    FILE_GZIPPED,
    FILE_MATE_1,
)


class TestSubsetFastq:
    """Test ``SubsetFastq`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        SubsetFastq(path=FILE_MATE_1)

    def test_init_all(self, tmpdir):
        """Create instance with all available parameters."""
        SubsetFastq(
            path=FILE_MATE_1,
            out_dir=tmpdir,
            records=1,
        )

    def test_process_no_process(self, tmpdir):
        """Input file is not processed."""
        test_instance = SubsetFastq(
            path=FILE_MATE_1,
            out_dir=tmpdir,
        )
        test_instance.process()
        assert test_instance.path != test_instance.out_path
        assert test_instance.n_processed != 0

    def test_process_gzipped_all_records(self, tmpdir):
        """Input file is gzipped and all records are to be processed."""
        test_instance = SubsetFastq(
            path=FILE_GZIPPED,
            out_dir=tmpdir,
        )
        test_instance.process()
        assert test_instance.path != test_instance.out_path
        assert test_instance.n_processed != 0

    def test_process_n_records(self, tmpdir):
        """Process exactly n records."""
        test_instance = SubsetFastq(
            path=FILE_MATE_1,
            out_dir=tmpdir,
            records=3,
        )
        test_instance.process()
        assert test_instance.path != test_instance.out_path
        assert test_instance.n_processed == 3

    def test_process_empty_file(self, tmpdir):
        """Process exactly n records."""
        test_instance = SubsetFastq(
            path=FILE_EMPTY,
            out_dir=tmpdir,
            records=3,
        )
        with pytest.raises(FileProblem):
            test_instance.process()
