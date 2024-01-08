"""Unit tests for module ``get_library_type.py``."""

import pytest

from htsinfer.exceptions import (
    FileProblem,
    InconsistentFastqIdentifiers,
    MetadataWarning,
    UnknownFastqIdentifier,
)
from htsinfer.get_library_type import (
    GetLibType,
    GetFastqType,
)
from htsinfer.models import (
    ResultsSource,
    ResultsType,
    Source,
    SeqIdFormats,
    StatesType,
    StatesTypeRelationship,
)
from tests.utils import (
    FILE_DUMMY,
    FILE_EMPTY,
    FILE_FASTA,
    FILE_INCONSISTENT_IDS_MIXED_UNKNOWN,
    FILE_INCONSISTENT_IDS_SINGLE_MATE,
    FILE_INCONSISTENT_IDS_SINGLE_OLD_NEW,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_IDS_NOT_MATCH_1,
    FILE_IDS_NOT_MATCH_2,
    FILE_TRANSCRIPTS,
    FILE_SINGLE,
    FILE_UNKNOWN_SEQ_ID,
    RaiseError,
    SEQ_ID_DUMMY,
    SEQ_ID_MATE_1,
    SEQ_ID_MATE_2,
    SEQ_ID_SINGLE,
    CONFIG,
    MAPPING,
)


class TestGetLibType:
    """Test ``GetLibType`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        CONFIG.args.path_2_processed = None
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        assert test_instance.path_1 == FILE_MATE_1

    def test_init_all(self):
        """Create instance with all available parameters."""
        CONFIG.args.path_2_processed = FILE_MATE_2
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        assert test_instance.path_1 == FILE_MATE_1
        assert test_instance.path_2 == FILE_MATE_2

    def test_evaluate_one_file(self):
        """Get library type for a single file."""
        CONFIG.args.path_2_processed = None
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        test_instance.evaluate()
        assert test_instance.results == ResultsType(
            file_1=StatesType.first_mate,
            file_2=StatesType.not_available,
            relationship=StatesTypeRelationship.not_available,
        )

    def test_evaluate_two_files(self):
        """Get library type for two files."""
        CONFIG.args.path_2_processed = FILE_MATE_2
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        test_instance.evaluate()
        assert test_instance.results == ResultsType(
            file_1=StatesType.first_mate,
            file_2=StatesType.second_mate,
            relationship=StatesTypeRelationship.split_mates,
        )

    def test_evaluate_mate_relationship_split_mates(self):
        """Test mate relationship evaluation logic with input files being
        mates of a paired-end library.
        """
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        test_instance.results.file_1 = StatesType.first_mate
        test_instance.results.file_2 = StatesType.second_mate
        test_instance._evaluate_mate_relationship(
            ids_1=["A", "B", "C"],
            ids_2=["A", "B", "C"],
        )
        assert (
            test_instance.results.relationship ==
            StatesTypeRelationship.split_mates
        )
        test_instance.results.file_1 = StatesType.second_mate
        test_instance.results.file_2 = StatesType.first_mate
        test_instance._evaluate_mate_relationship(
            ids_1=["A", "B", "C"],
            ids_2=["A", "B", "C"],
        )
        assert (
            test_instance.results.relationship ==
            StatesTypeRelationship.split_mates
        )

    def test_evaluate_mate_relationship_not_mates(self, tmpdir):
        """Test mate relationship evaluation logic with input files that are
        mates, but the relationship is not enough to trigger split_mates.
        """
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.paths = (FILE_MATE_1, FILE_MATE_2)
        MAPPING.transcripts_file = FILE_TRANSCRIPTS
        MAPPING.tmp_dir = tmpdir

        test_instance = GetLibType(config=CONFIG, mapping=MAPPING)
        test_instance.results.file_1 = StatesType.not_available
        test_instance.results.file_2 = StatesType.not_available

        # Set the cutoff such that it's not enough to trigger split_mates
        test_instance.cutoff = 300

        # Call the _evaluate_mate_relationship method
        test_instance._evaluate_mate_relationship(
            ids_1=["A", "B", "C"], ids_2=["A", "B", "C"]
        )

        assert (
            test_instance.results.relationship ==
            StatesTypeRelationship.not_mates
        )

    def test_evaluate_mate_relationship_not_available(self, tmpdir):
        """Test mate relationship evaluation logic with input files that are
        not mates from a paired-end library.
        """
        CONFIG.args.path_1_processed = FILE_IDS_NOT_MATCH_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.results.library_source = ResultsSource(
            file_1=Source(short_name="hsapiens", taxon_id=9606),
            file_2=Source(short_name="hsapiens", taxon_id=9606),
        )
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.paths = (FILE_IDS_NOT_MATCH_1, FILE_MATE_2)
        MAPPING.transcripts_file = FILE_TRANSCRIPTS
        MAPPING.tmp_dir = tmpdir
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        test_instance.results.file_1 = StatesType.not_available
        test_instance.results.file_2 = StatesType.not_available
        test_instance.evaluate()
        assert (
            test_instance.results.relationship ==
            StatesTypeRelationship.not_available
        )

    def test_evaluate_split_mates_not_matching_ids(self, tmpdir):
        """Test mate relationship evaluation logic with input files that are
        not mates from a paired-end library.
        """
        CONFIG.args.path_1_processed = FILE_IDS_NOT_MATCH_1
        CONFIG.args.path_2_processed = FILE_IDS_NOT_MATCH_2
        CONFIG.results.library_source = ResultsSource(
            file_1=Source(short_name="hsapiens", taxon_id=9606),
            file_2=Source(short_name="hsapiens", taxon_id=9606),
        )
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.paths = (FILE_IDS_NOT_MATCH_1, FILE_IDS_NOT_MATCH_2)
        MAPPING.tmp_dir = tmpdir
        test_instance = GetLibType(config=CONFIG,
                                   mapping=MAPPING)
        test_instance.evaluate()
        assert (
                test_instance.results.relationship ==
                StatesTypeRelationship.split_mates
        )


class TestAlignedSegment:
    """Test ``GetFastqType`` class."""

    def test_init(self):
        """Create instance."""
        test_instance = GetLibType.AlignedSegment()
        assert test_instance.reference_name == ""
        assert test_instance.reference_start == 0
        assert test_instance.reference_end == 0


class TestGetFastqType:
    """Test ``GetFastqType`` class."""

    def test_init(self):
        """Create instance."""
        test_instance = GetFastqType(path=FILE_MATE_1)
        assert test_instance.path == FILE_MATE_1

    def test_evaluate_single(self):
        """Evaluate valid single-end library file."""
        test_instance = GetFastqType(path=FILE_SINGLE)
        test_instance.evaluate()
        assert test_instance.result == StatesType.single

    def test_evaluate_inconsistent_identifiers_single_mate(self):
        """Raise ``MetadataWarning`` by passing a file with inconsistent
        identifiers, suggesting a single-end library first, then a paired-end
        library later.
        """
        test_instance = GetFastqType(path=FILE_INCONSISTENT_IDS_SINGLE_MATE)
        with pytest.raises(MetadataWarning):
            test_instance.evaluate()

    def test_evaluate_inconsistent_identifiers_mixed_unknown(self):
        """Raise ``MetadataWarning`` by passing a file with inconsistent
        identifiers, suggesting a mixed paired-end library first, followed by
        reads of an unknown identifier format.
        """
        test_instance = GetFastqType(path=FILE_INCONSISTENT_IDS_MIXED_UNKNOWN)
        with pytest.raises(MetadataWarning):
            test_instance.evaluate()

    def test_evaluate_inconsistent_identifiers_single_old_new(self):
        """Raise ``MetadataWarning`` by passing a file with inconsistent
        identifier formats (Cavatica <1.8 vs Cavatica >=1.8).
        """
        test_instance = GetFastqType(path=FILE_INCONSISTENT_IDS_SINGLE_OLD_NEW)
        with pytest.raises(MetadataWarning):
            test_instance.evaluate()

    def test_evaluate_file_problem_empty_file(self):
        """Pass empty file to simulate a file problem."""
        test_instance = GetFastqType(path=FILE_EMPTY)
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_fasta(self):
        """Pass FASTA file to FASTQ parser (raising a ``ValueError``) to
        simulate a file problem.
        """
        test_instance = GetFastqType(path=FILE_FASTA)
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_evaluate_file_problem_cannot_open_file(self, monkeypatch):
        """Force raising of ``OSError`` to simulate file problem."""
        test_instance = GetFastqType(path=FILE_SINGLE)
        monkeypatch.setattr(
            'builtins.open',
            RaiseError(exc=OSError),
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_get_read_type_match_single(self):
        """Evaluate read consistent with single-end library."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance._get_read_type(
            seq_id=SEQ_ID_SINGLE,
            regex=SeqIdFormats['Casava <1.8'].value,
        )
        assert test_instance.result == StatesType.single

    def test_get_read_type_match_first_mate(self):
        """Evaluate read consistent with first mate file of a paired-end
        library.
        """
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance._get_read_type(
            seq_id=SEQ_ID_MATE_1,
            regex=SeqIdFormats['Casava <1.8'].value,
        )
        assert test_instance.result == StatesType.first_mate

    def test_get_read_type_match_second_mate(self):
        """Evaluate read consistent with second mate file of a paired-end
        library.
        """
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance._get_read_type(
            seq_id=SEQ_ID_MATE_2,
            regex=SeqIdFormats['Casava <1.8'].value,
        )
        assert test_instance.result == StatesType.second_mate

    def test_get_read_type_no_match(self):
        """Evaluate read that does not match the indicated regular expression.
        """
        test_instance = GetFastqType(path=FILE_DUMMY)
        with pytest.raises(UnknownFastqIdentifier):
            test_instance._get_read_type(
                seq_id=SEQ_ID_SINGLE,
                regex=SeqIdFormats['Casava >=1.8'].value,
            )

    def test_get_read_type_single_pass(self):
        """Read identifier is consistent with previous state."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.single
        test_instance._get_read_type_single(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.single

    def test_get_read_type_single_set(self):
        """Set library type based on read identifier for the first time."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.not_available
        test_instance._get_read_type_single(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.single

    def test_get_read_type_single_inconsistent(self):
        """Read identifier is not consistent with previous state."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.first_mate
        with pytest.raises(InconsistentFastqIdentifiers):
            test_instance._get_read_type_single(seq_id=SEQ_ID_DUMMY)

    def test_get_read_type_mate_1_pass(self):
        """Read identifier is consistent with previous state."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.mixed_mates
        test_instance._get_read_type_paired_mate_1(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.mixed_mates

    def test_get_read_type_mate_1_mixed(self):
        """Read identifier is consistent with previous state, but suggests
        mixed paired-end library file.
        """
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.second_mate
        test_instance._get_read_type_paired_mate_1(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.mixed_mates

    def test_get_read_type_mate_1_set(self):
        """Set library type based on read identifier for the first time."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.not_available
        test_instance._get_read_type_paired_mate_1(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.first_mate

    def test_get_read_type_mate_1_inconsistent(self):
        """Read identifier is not consistent with previous state."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.single
        with pytest.raises(InconsistentFastqIdentifiers):
            test_instance._get_read_type_paired_mate_1(seq_id=SEQ_ID_DUMMY)

    def test_get_read_type_mate_2_pass(self):
        """Read identifier is consistent with previous state."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.mixed_mates
        test_instance._get_read_type_paired_mate_2(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.mixed_mates

    def test_get_read_type_mate_2_mixed(self):
        """Read identifier is consistent with previous state, but suggests
        mixed paired-end library file.
        """
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.first_mate
        test_instance._get_read_type_paired_mate_2(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.mixed_mates

    def test_get_read_type_mate_2_set(self):
        """Set library type based on read identifier for the first time."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.not_available
        test_instance._get_read_type_paired_mate_2(seq_id=SEQ_ID_DUMMY)
        assert test_instance.result == StatesType.second_mate

    def test_get_read_type_mate_2_inconsistent(self):
        """Read identifier is not consistent with previous state."""
        test_instance = GetFastqType(path=FILE_DUMMY)
        test_instance.result = StatesType.single
        with pytest.raises(InconsistentFastqIdentifiers):
            test_instance._get_read_type_paired_mate_2(seq_id=SEQ_ID_DUMMY)
