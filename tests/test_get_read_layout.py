"""Unit tests for module ``get_read_layout.py``."""

from pathlib import Path
import pytest

from htsinfer.exceptions import FileProblem
from htsinfer.get_read_layout import (
    GetReadLayout,
    GetAdapter3,
)
from htsinfer.models import (
    Layout,
    ResultsLayout,
)
from tests.utils import (
    FILE_ADAPTER,
    FILE_DUMMY,
    FILE_FASTA,
    FILE_INVALID_SEQ_1,
    FILE_INVALID_SEQ_2,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_REAL_SAMPLE,
    FILE_SINGLE,
    FILE_SRA_SAMPLE_1,
    FILE_SRA_SAMPLE_2,
    RaiseOSError,
)


class TestGetReadLayout:
    """Test ``GetReadLayout`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        test_instance = GetReadLayout(path_1=FILE_MATE_1)
        assert test_instance.path_1 == FILE_MATE_1

    def test_init_all(self):
        """Create instance with all available parameters."""
        test_instance = GetReadLayout(
            path_1=FILE_MATE_1,
            path_2=FILE_MATE_2,
            adapter_file=FILE_ADAPTER,
            out_dir=Path.cwd(),
            min_match_pct=5,
            min_freq_ratio=2,
        )
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.path_2 == FILE_MATE_2
        assert test_instance.adapter_file == FILE_ADAPTER
        assert test_instance.out_dir == Path.cwd()
        assert test_instance.min_match_pct == 5
        assert test_instance.min_freq_ratio == 2

    def test_evaluate_one_file(self):
        """Get read layout for a single file."""
        test_instance = GetReadLayout(
            path_1=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER,
        )
        test_instance.evaluate()
        assert test_instance.results == ResultsLayout(
            file_1=Layout(adapt_3="AAAAAAAAAAAAAAA"),
            file_2=Layout(adapt_3=None),
        )

    def test_evaluate_two_files(self):
        """Get read layout for two files."""
        test_instance = GetReadLayout(
            path_1=FILE_SRA_SAMPLE_1,
            path_2=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER,
            min_match_pct=2,
            min_freq_ratio=1,
        )
        test_instance.evaluate()
        assert test_instance.results == ResultsLayout(
            file_1=Layout(adapt_3="AAAAAAAAAAAAAAA"),
            file_2=Layout(adapt_3="AAAAAAAAAAAAAAA"),
        )


class TestGetAdapter3:
    """Test ``GetAdapter3`` class."""

    def test_init(self):
        """Create instance."""
        test_instance = GetAdapter3(
            path=FILE_MATE_1,
            adapter_file=FILE_ADAPTER,
            out_dir=Path.cwd(),
            min_match_pct=5,
            min_freq_ratio=2,
        )
        assert test_instance.path == FILE_MATE_1
        assert test_instance.adapter_file == FILE_ADAPTER
        assert test_instance.out_dir == Path.cwd()
        assert test_instance.min_match_pct == 5
        assert test_instance.min_freq_ratio == 2

    def test_evaluate_layout(self):
        """Evaluate valid read layout file."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER,
        )
        test_instance.evaluate()
        assert test_instance.result == "AAAAAAAAAAAAAAA"

    def test_evaluate_dummy_file(self):
        """Pass dummy file to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_DUMMY,
            adapter_file=FILE_ADAPTER,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evualte_invalid_sequence_1(self):
        """Pass a file with invalid sequence to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_INVALID_SEQ_1,
            adapter_file=FILE_ADAPTER,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_cannot_open_file(self, monkeypatch):
        """Force raising of ``OSError`` to simulate file problem."""
        test_instance = GetAdapter3(
            path=FILE_SINGLE,
            adapter_file=FILE_ADAPTER,
        )
        monkeypatch.setattr(
            'builtins.open',
            RaiseOSError,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_fasta(self):
        """Pass FASTA file to FASTQ parser (raising a ``ValueError``) to
        simulate a file problem.
        """
        test_instance = GetAdapter3(
            path=FILE_FASTA,
            adapter_file=FILE_ADAPTER,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_adapter(self):
        """Pass dummy file as adapter file to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_SINGLE,
            adapter_file=FILE_DUMMY,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_results_validator(self):
        """Pass valid file with default values of min_match_pct and
        min_freq_ratio to evaluate results validator."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_1,
            adapter_file=FILE_ADAPTER,
        )
        test_instance.evaluate()
        assert test_instance.result is None

    def test_evaluate_min_freq_ratio(self):
        """Pass valid file with different min_freq_ratio value to evaluate
        results validator."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER,
            min_freq_ratio=4,
        )
        test_instance.evaluate()
        assert test_instance.result is None

    def test_evualte_invalid_sequence_2(self):
        """Pass a file with invalid sequence to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_INVALID_SEQ_2,
            adapter_file=FILE_ADAPTER,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_no_read_layout(self):
        """Pass a file to test if no read layout found."""
        test_instance = GetAdapter3(
            path=FILE_REAL_SAMPLE,
            adapter_file=FILE_ADAPTER,
        )
        test_instance.evaluate()
        assert test_instance.result is None
