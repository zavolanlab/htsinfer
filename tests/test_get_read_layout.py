"""Unit tests for module ``get_read_layout.py``."""

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
    FILE_ADAPTER_INVALID_CHARS,
    FILE_ADAPTER_SEQ_TOO_LONG,
    FILE_DUMMY,
    FILE_FASTA,
    FILE_INVALID_SEQ_1,
    FILE_INVALID_SEQ_2,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_2000_RECORDS,
    FILE_SINGLE,
    FILE_SRA_SAMPLE_1,
    FILE_SRA_SAMPLE_2,
    RaiseError,
    CONFIG,
)


class TestGetReadLayout:
    """Test ``GetReadLayout`` class."""

    def test_init_required(self, tmpdir):
        """Create instance with required parameters."""
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.path_2_processed = None
        test_instance = GetReadLayout(config=CONFIG)
        assert test_instance.path_1 == FILE_MATE_1

    def test_init_all(self, tmpdir):
        """Create instance with all available parameters."""
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.read_layout_adapter_file = FILE_ADAPTER
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.read_layout_min_match_pct = 5
        CONFIG.args.read_layout_min_freq_ratio = 2
        test_instance = GetReadLayout(config=CONFIG)
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.path_2 == FILE_MATE_2
        assert test_instance.adapter_file == FILE_ADAPTER
        assert test_instance.out_dir == tmpdir
        assert test_instance.min_match_pct == 5
        assert test_instance.min_freq_ratio == 2

    def test_evaluate_one_file(self, tmpdir):
        """Get read layout for a single file."""
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.path_1_processed = FILE_SRA_SAMPLE_2
        CONFIG.args.path_2_processed = None
        test_instance = GetReadLayout(config=CONFIG)
        test_instance.evaluate()
        assert test_instance.results == ResultsLayout(
            file_1=Layout(adapt_3="AAAAAAAAAAAAAAA"),
            file_2=Layout(adapt_3=None),
        )

    def test_evaluate_two_files(self, tmpdir):
        """Get read layout for two files."""
        CONFIG.args.path_1_processed = FILE_SRA_SAMPLE_2
        CONFIG.args.path_2_processed = FILE_SRA_SAMPLE_2
        CONFIG.args.read_layout_adapter_file = FILE_ADAPTER
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.lib_source_min_match_pct = 2
        CONFIG.args.read_layout_min_freq_ratio = 1
        test_instance = GetReadLayout(config=CONFIG)
        test_instance.evaluate()
        assert test_instance.results == ResultsLayout(
            file_1=Layout(adapt_3="AAAAAAAAAAAAAAA"),
            file_2=Layout(adapt_3="AAAAAAAAAAAAAAA"),
        )


class TestGetAdapter3:
    """Test ``GetAdapter3`` class."""

    def test_init(self, tmpdir):
        """Create instance."""
        tmp_dir = tmpdir
        test_instance = GetAdapter3(
            path=FILE_MATE_1,
            adapter_file=FILE_ADAPTER,
            out_dir=tmp_dir,
            min_match_pct=5,
            min_freq_ratio=2,
        )
        assert test_instance.path == FILE_MATE_1
        assert test_instance.adapter_file == FILE_ADAPTER
        assert test_instance.out_dir == tmp_dir
        assert test_instance.min_match_pct == 5
        assert test_instance.min_freq_ratio == 2

    def test_evaluate_layout(self, tmpdir):
        """Evaluate valid read layout file."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.result == "AAAAAAAAAAAAAAA"

    def test_evaluate_dummy_file(self, tmpdir):
        """Pass dummy file to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_DUMMY,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evualte_invalid_sequence_1(self, tmpdir):
        """Pass a file with invalid sequence to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_INVALID_SEQ_1,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_cannot_open_file(self, monkeypatch, tmpdir):
        """Force raising of ``OSError`` to simulate file problem."""
        test_instance = GetAdapter3(
            path=FILE_SINGLE,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'builtins.open',
            RaiseError(exc=OSError),
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_fasta(self, tmpdir):
        """Pass FASTA file to FASTQ parser (raising a ``ValueError``) to
        simulate a file problem.
        """
        test_instance = GetAdapter3(
            path=FILE_FASTA,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_adapter(self, tmpdir):
        """Pass dummy file as adapter file to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_SINGLE,
            adapter_file=FILE_DUMMY,
            out_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_results_validator(self, tmpdir):
        """Pass valid file with default values of min_match_pct and
        min_freq_ratio to evaluate results validator."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_1,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.result is None

    def test_evaluate_min_freq_ratio(self, tmpdir):
        """Pass valid file with different min_freq_ratio value to evaluate
        results validator."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER,
            min_freq_ratio=4,
            out_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.result is None

    def test_evaluate_invalid_sequence_2(self, tmpdir):
        """Pass a file with invalid sequence to simulate a file problem."""
        test_instance = GetAdapter3(
            path=FILE_INVALID_SEQ_2,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_no_read_layout(self, tmpdir):
        """Pass a file to test if no read layout found."""
        test_instance = GetAdapter3(
            path=FILE_2000_RECORDS,
            adapter_file=FILE_ADAPTER,
            out_dir=tmpdir,
        )
        test_instance.evaluate()
        assert test_instance.result is None

    def test_evaluate_invalid_adapter_sequence_chars(self, tmpdir):
        """Pass an adapter file with sequence containing illegal characters."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER_INVALID_CHARS,
            out_dir=tmpdir,
        )
        with pytest.raises(ValueError):
            test_instance.evaluate()

    def test_evaluate_invalid_adapter_sequence_too_long(self, tmpdir):
        """Pass an adapter file with sequence that are too long."""
        test_instance = GetAdapter3(
            path=FILE_SRA_SAMPLE_2,
            adapter_file=FILE_ADAPTER_SEQ_TOO_LONG,
            out_dir=tmpdir,
        )
        with pytest.raises(ValueError):
            test_instance.evaluate()
