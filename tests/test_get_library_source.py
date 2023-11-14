"""Unit tests for module ``get_read_library_source.py``."""

import pytest

from htsinfer.exceptions import (
    FileProblem,
    KallistoProblem,
    TranscriptsFastaProblem
)
from htsinfer.get_library_source import GetLibSource
from htsinfer.models import (
    ResultsSource,
    Source,
)
from tests.utils import (
    FILE_2000_RECORDS,
    FILE_DUMMY,
    FILE_EMPTY,
    FILE_INVALID_PATH,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_SOURCE_FRUIT_FLY,
    FILE_TRANSCRIPTS,
    SOURCE_FRUIT_FLY,
    CONFIG,
    SOURCE_HUMAN,
    TEST_FILES_DIR,
)


class TestGetLibSource:
    """Test ``GetLibSource`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        CONFIG.args.path_2_processed = None
        test_instance = GetLibSource(
            config=CONFIG,
        )
        assert test_instance.paths[0] == FILE_MATE_1

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        CONFIG.args.path_2_processed = FILE_MATE_2
        test_instance = GetLibSource(
            config=CONFIG,
        )
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2

    def test_evaluate_single(self, monkeypatch, tmpdir):
        """Get library statistics for a single-end library."""
        CONFIG.args.path_2_processed = None
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        test_instance = GetLibSource(
            config=CONFIG,
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
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        CONFIG.args.path_2_processed = FILE_MATE_2
        test_instance = GetLibSource(
            config=CONFIG,
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
        CONFIG.args.path_1_processed = FILE_2000_RECORDS
        CONFIG.args.path_2_processed = None
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(
            config=CONFIG,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=SOURCE_HUMAN,
            file_2=Source()
        )

    def test_evaluate_paired_different_source(self, tmpdir):
        """Pass both files with source fruit fly and human."""
        CONFIG.args.path_1_processed = FILE_SOURCE_FRUIT_FLY
        CONFIG.args.path_2_processed = FILE_2000_RECORDS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(
            config=CONFIG,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=SOURCE_FRUIT_FLY,
            file_2=SOURCE_HUMAN
        )

    def test_evaluate_no_library_source(self, tmpdir):
        """Pass a file to test if no library source found."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = None
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(
            config=CONFIG,
        )
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_paired_no_library_source(self, tmpdir):
        """Pass paired files to test if library source found."""
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_dummy_file(self, tmpdir):
        """Pass dummy file to test."""
        CONFIG.args.path_1_processed = FILE_DUMMY
        CONFIG.args.path_2_processed = None
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_empty_file(self, tmpdir):
        """Pass empty file to test."""
        CONFIG.args.path_1_processed = FILE_EMPTY
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_file_problem_transcripts(self, tmpdir):
        """Pass dummy file as transcripts.fasta file
        to simulate value error."""
        CONFIG.args.path_1_processed = FILE_2000_RECORDS
        CONFIG.args.t_file_processed = FILE_DUMMY
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        with pytest.raises(TranscriptsFastaProblem):
            test_instance.evaluate()

    def test_evaluate_kallisto_index_problem(self, monkeypatch, tmpdir):
        """Force raising exception to simulate KallistoProblem."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.create_kallisto_index',
            KallistoProblem
            )
        with pytest.raises(KallistoProblem):
            test_instance.evaluate()

    def test_evaluate_kallisto_quant_problem(self, monkeypatch, tmpdir):
        """Force raising exception to stimulate KallistoProblem."""
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        sub_method_name = 'htsinfer.get_library_source.' + \
            'GetLibSource.run_kallisto_quantification'
        monkeypatch.setattr(
            sub_method_name,
            lambda *args, **kwargs: KallistoProblem,
            )
        with pytest.raises(TypeError):
            test_instance.evaluate()

    def test_evaluate_get_source_expression_problem(self, monkeypatch, tmpdir):
        """Force raising exception to stimulate a file probblem."""
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        sub_method_name = 'htsinfer.get_library_source.' + \
            'GetLibSource.run_kallisto_quantification'
        monkeypatch.setattr(
            sub_method_name,
            lambda *args, **kwargs: tmpdir,
            )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_get_source_expression_abundace_file_problem(
        self, monkeypatch, tmpdir
    ):
        """Pass empty abundance.tsv file to simulate value error."""
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        sub_method_name = 'htsinfer.get_library_source.' + \
            'GetLibSource.run_kallisto_quantification'
        monkeypatch.setattr(
            sub_method_name,
            lambda *args, **kwargs: TEST_FILES_DIR,
            )
        with pytest.raises(TranscriptsFastaProblem):
            test_instance.evaluate()

    def test_evaluate_min_match_pct(self, tmpdir):
        """Pass file with minimum match percentage to
        test validate_top_score."""
        CONFIG.args.lib_source_min_match_pct = 99
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_min_freq_ratio(self, tmpdir):
        """Pass file with minimum frequency ratio to
        test validate_top_score."""
        CONFIG.args.lib_source_min_match_pct = 2
        CONFIG.args.lib_source_min_freq_ratio = 15
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        results = test_instance.evaluate()
        assert results == ResultsSource(
            file_1=Source(),
            file_2=Source()
        )

    def test_evaluate_tax_id_not_none(self):
        """Test when self.tax_id is not None."""
        CONFIG.args.tax_id = 7227  # An example taxon ID
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        test_instance = GetLibSource(config=CONFIG)
        result = test_instance.evaluate()

        assert result.file_1.taxon_id == 7227
        assert result.file_1.short_name == "dmelanogaster"

    def test_evaluate_tax_id_none_with_path_2(self, tmpdir, monkeypatch):
        """Test when self.tax_id is None and self.paths[1] is not None."""
        CONFIG.args.tax_id = None
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)

        # Mock the get_source method to return a specific result
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_source',
            lambda *args, **kwargs: SOURCE_HUMAN,
        )

        result = test_instance.evaluate()

        assert result.file_2.taxon_id == SOURCE_HUMAN.taxon_id
        assert result.file_2.short_name == SOURCE_HUMAN.short_name

    def test_evaluate_tax_id_not_none_with_path_2(self, tmpdir):
        """Test when self.tax_id is not None and self.paths[1] is not None."""
        CONFIG.args.tax_id = 7227
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)

        result = test_instance.evaluate()

        assert result.file_2.taxon_id == 7227
        assert result.file_2.short_name == "dmelanogaster"

    def test_create_kallisto_index_problem(self, tmpdir):
        """Pass invalid file as transcripts.fasta file
        to simulate KallistoProblem."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.lib_source_min_freq_ratio = 2
        CONFIG.args.t_file_processed = FILE_INVALID_PATH
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)
        with pytest.raises(KallistoProblem):
            test_instance.create_kallisto_index()

    def test_get_organism_name_found(self):
        """Test the function when the taxon_id
        is found in the organism dictionary."""
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        test_instance = GetLibSource(config=CONFIG)
        taxon_id = 7227
        result = test_instance.get_organism_name(
            taxon_id, CONFIG.args.t_file_processed
        )
        assert result == "dmelanogaster"

    def test_get_organism_name_not_found(self):
        """Test the function when the taxon_id
        is not found in the organism dictionary."""
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        test_instance = GetLibSource(config=CONFIG)
        taxon_id = 12345  # A tax ID that doesn't exist in transcripts
        result = test_instance.get_organism_name(
            taxon_id, CONFIG.args.t_file_processed
        )
        assert result is None

    def test_get_organism_name_file_problem(self):
        """Test the function when there's a
        file problem while processing the FASTA file."""
        CONFIG.args.t_file_processed = FILE_DUMMY
        test_instance = GetLibSource(config=CONFIG)
        taxon_id = 7227
        with pytest.raises(FileProblem):
            test_instance.get_organism_name(
                taxon_id, CONFIG.args.t_file_processed
            )

    def test_evaluate_tax_id_is_none(self, monkeypatch, tmpdir):
        """Test when self.tax_id is None."""
        CONFIG.args.tax_id = None
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)

        # Mock the create_kallisto_index method to return a specific result
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.create_kallisto_index',
            lambda *args, **kwargs: tmpdir / "kallisto.idx",
        )

        # Mock the get_source method to return a specific result
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_source',
            lambda *args, **kwargs: SOURCE_FRUIT_FLY,
        )

        result = test_instance.evaluate()

        assert result.file_1.taxon_id == SOURCE_FRUIT_FLY.taxon_id
        assert result.file_1.short_name == SOURCE_FRUIT_FLY.short_name

        assert result.file_2.taxon_id == SOURCE_FRUIT_FLY.taxon_id
        assert result.file_2.short_name == SOURCE_FRUIT_FLY.short_name

    def test_evaluate_tax_id_not_none_no_src_name(self, monkeypatch, tmpdir):
        """Test when self.tax_id is not None but src_name is not found."""
        CONFIG.args.tax_id = 7227
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)

        # Mock the get_organism_name method to return None
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_organism_name',
            lambda *args, **kwargs: None,
        )

        # Mock the create_kallisto_index method to return a specific result
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.create_kallisto_index',
            lambda *args, **kwargs: tmpdir / "kallisto.idx",
        )

        # Mock the get_source method to return a specific result
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_source',
            lambda *args, **kwargs: SOURCE_FRUIT_FLY,
        )

        result = test_instance.evaluate()

        assert result.file_1.taxon_id == SOURCE_FRUIT_FLY.taxon_id
        assert result.file_1.short_name == SOURCE_FRUIT_FLY.short_name

        assert result.file_2.taxon_id == SOURCE_FRUIT_FLY.taxon_id
        assert result.file_2.short_name == SOURCE_FRUIT_FLY.short_name

    def test_evaluate_tax_id_not_none_name_found(self, monkeypatch, tmpdir):
        """Test when self.tax_id is not None and src_name is found."""
        CONFIG.args.tax_id = 7227
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.args.out_dir = tmpdir
        test_instance = GetLibSource(config=CONFIG)

        # Mock the get_organism_name method to return a specific result
        monkeypatch.setattr(
            'htsinfer.get_library_source.GetLibSource.get_organism_name',
            lambda *args, **kwargs: "dmelanogaster",
        )

        result = test_instance.evaluate()

        assert result.file_1.taxon_id == 7227
        assert result.file_1.short_name == "dmelanogaster"

        assert result.file_2.taxon_id == 7227
        assert result.file_2.short_name == "dmelanogaster"
