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
    ADAPTER_FILE,
    FILE_DUMMY,
    FILE_FASTA,
    FILE_INVALID_SEQ,
    FILE_MATE_1,
    FILE_MATE_2,
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
            adapter_file=ADAPTER_FILE,
            out_dir=Path.cwd(),
            min_match_pct=5,
            min_freq_ratio=2,
        )
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.path_2 == FILE_MATE_2
        assert test_instance.adapter_file == ADAPTER_FILE
        assert test_instance.out_dir == Path.cwd()
        assert test_instance.min_match_pct == 5
        assert test_instance.min_freq_ratio == 2

    def test_evaluate_one_file(self):
        """Get read layout for a single file."""
        test_instance = GetReadLayout(
            path_1=FILE_SRA_SAMPLE_2,
            adapter_file=ADAPTER_FILE,
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
            adapter_file=ADAPTER_FILE,
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
            adapter_file=ADAPTER_FILE,
            out_dir=Path.cwd(),
            min_match_pct=5,
            min_freq_ratio=2,
        )
        assert test_instance.path == FILE_MATE_1
        assert test_instance.adapter_file == ADAPTER_FILE
        assert test_instance.out_dir == Path.cwd()
        assert test_instance.min_match_pct == 5
        assert test_instance.min_freq_ratio == 2

    def test_evaluate_layout(self):
        """Evaluate valid read layout file."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=ADAPTER_FILE,
        )
        test_instance.evaluate()
        assert test_instance.result == "AAAAAAAAAAAAAAA"

    def test_evaluate_dummy_file(self):
        """Pass dummy file to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_DUMMY,
            adapter_file=ADAPTER_FILE,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evualte_invalid_sequence(self):
        """Pass a file with invalid sequence to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_INVALID_SEQ,
            adapter_file=ADAPTER_FILE,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_cannot_open_file(self, monkeypatch):
        """Force raising of ``OSError`` to simulate file problem."""
        test_instance = GetAdapter3(
            path=FILE_SINGLE,
            adapter_file=ADAPTER_FILE,
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
            adapter_file=ADAPTER_FILE,
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
