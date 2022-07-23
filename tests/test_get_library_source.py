"""Unit tests for module ``get_read_library_source.py``."""

import pytest

from htsinfer.exceptions import FileProblem, KallistoProblem
from htsinfer.get_library_source import GetLibSource
from htsinfer.models import (
    ResultsSource,
    Source,
)
from tests.utils import (
    FILE_2000_RECORDS,
    FILE_DUMMY,
    FILE_EMPTY,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_SOURCE_FRUIT_FLY,
    FILE_TRANSCRIPTS,
    SOURCE_FRUIT_FLY,
    SOURCE_HUMAN,
    TEST_FILES_DIR,
)


class TestGetLibSource:
    """Test ``GetLibSource`` class."""

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

    def test_evaluate_source_human(self, tmpdir):
        """Pass file with source human."""
        test_instance = GetLibSource(
            paths=(FILE_2000_RECORDS, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=SOURCE_HUMAN,
            file_2=Source()
        )

    def test_evaluate_paired_different_source(self, tmpdir):
        """Pass both files with source fruit fly and human."""
        test_instance = GetLibSource(
            paths=(FILE_SOURCE_FRUIT_FLY, FILE_2000_RECORDS),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=SOURCE_FRUIT_FLY,
            file_2=SOURCE_HUMAN
        )

    def test_evaluate_no_library_source(self, tmpdir):
        """Pass a file to test if no library source found."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_paired_no_library_source(self, tmpdir):
        """Pass paired files to test if library source found."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, FILE_MATE_2),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_dummy_file(self, tmpdir):
        """Pass dummy file to test."""
        test_instance = GetLibSource(
            paths=(FILE_DUMMY, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_empty_file(self, tmpdir):
        """Pass empty file to test."""
        test_instance = GetLibSource(
            paths=(FILE_EMPTY, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evualte_file_problem_transcripts(self, tmpdir):
        """Pass dummy file as transcripts.fasta file
        to simulate a file problem."""
        test_instance = GetLibSource(
            paths=(FILE_2000_RECORDS, None),
            transcripts_file=FILE_DUMMY,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_kallisto_index_problem(self, monkeypatch, tmpdir):
        """Force raising exception to simulate KallistoProblem."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.\
                GetLibSource.create_kallisto_index',
            KallistoProblem
            )
        with pytest.raises(KallistoProblem):
            test_instance.evaluate()

    def test_evaluate_kallisto_quant_problem(self, monkeypatch, tmpdir):
        """Force raising exception to stimulate KallistoProblem."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.\
                GetLibSource.run_kallisto_quantification',
            lambda *args, **kwargs: KallistoProblem,
            )
        with pytest.raises(TypeError):
            test_instance.evaluate()

    def test_evaluate_get_source_expression_problem(self, monkeypatch, tmpdir):
        """Force raising exception to stimulate a file probblem."""
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.\
                GetLibSource.run_kallisto_quantification',
            lambda *args, **kwargs: tmpdir,
            )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_get_source_expression_abundace_file_problem(
        self, monkeypatch, tmpdir
    ):
        """Pass empty abundance.tsv file (raising AttributeError)
        to simulate a file problem.
        """
        test_instance = GetLibSource(
            paths=(FILE_MATE_1, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_library_source.\
                GetLibSource.run_kallisto_quantification',
            lambda *args, **kwargs: TEST_FILES_DIR,
            )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_min_match_pct(self, tmpdir):
        """Pass file with minimum match percentage to
        test validate_top_score."""
        test_instance = GetLibSource(
            paths=(FILE_SOURCE_FRUIT_FLY, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
            min_match_pct=99,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_min_freq_ratio(self, tmpdir):
        """Pass file with minimum frequency ratio to
        test validate_top_score."""
        test_instance = GetLibSource(
            paths=(FILE_2000_RECORDS, None),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
            out_dir=tmpdir,
            min_freq_ratio=5,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )
