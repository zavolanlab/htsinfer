"""Unit tests for module ``get_read_orientation.py``."""

import pytest

from htsinfer.exceptions import (
    FileProblem,
)
from htsinfer.get_read_orientation import GetOrientation
from htsinfer.models import (
    ResultsOrientation,
    ResultsSource,
    ResultsType,
    Source,
    StatesOrientation,
    StatesOrientationRelationship,
    StatesTypeRelationship,
)
from tests.utils import (
    FILE_EMPTY_ALIGNED_SAM,
    FILE_BAD_ALIGNED_SAM,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_ORIENTATION_ISF_1,
    FILE_ORIENTATION_ISF_2,
    FILE_ORIENTATION_ISR_1,
    FILE_ORIENTATION_ISR_2,
    FILE_ORIENTATION_IU_1,
    FILE_ORIENTATION_IU_2,
    FILE_ORIENTATION_SF,
    FILE_ORIENTATION_SR,
    FILE_ORIENTATION_U,
    FILE_TRANSCRIPTS,
    FILE_UNMAPPED_PAIRED_1,
    FILE_UNMAPPED_PAIRED_2,
    FILE_UNMAPPED_SINGLE,
    FILE_IDS_NOT_MATCH_2,
    CONFIG,
    MAPPING,
)


class TestGetOrientation:
    """Test ``GetOrientation`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = None
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.results.library_type = ResultsType()
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.library_type == ResultsType()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2
        assert test_instance.library_type == ResultsType()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS

    def test_init_all(self, tmpdir):
        """Create instance with all available parameters."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.results.library_source = ResultsSource(
            file_1=Source(),
            file_2=Source()
        )
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2
        assert test_instance.library_type == ResultsType()
        assert test_instance.library_source == ResultsSource()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS
        assert test_instance.tmp_dir == tmpdir
        assert test_instance.min_mapped_reads == 18
        assert test_instance.min_fraction == 0.75

    def test_evaluate_single_unmapped(self, tmpdir):
        """Get read orientation for a single-end library with no mappable
        reads.
        """
        CONFIG.args.path_1_processed = FILE_UNMAPPED_SINGLE
        CONFIG.args.path_2_processed = None
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.not_available,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_single_sf(self, tmpdir):
        """Get read orientation for a single-end stranded forward library."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_SF
        CONFIG.results.library_source = ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(),
            )
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_forward,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_single_sr(self, tmpdir):
        """Get read orientation for a single-end stranded reverse library."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_SR
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_reverse,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_single_u(self, tmpdir):
        """Get read orientation for a single-end unstranded library."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_U
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.unstranded,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_paired_unmapped(self, tmpdir):
        """Get read orientation for a paired-end library with no mappable
        reads.
        """
        CONFIG.args.path_1_processed = FILE_UNMAPPED_PAIRED_1
        CONFIG.args.path_2_processed = FILE_UNMAPPED_PAIRED_2
        CONFIG.args.tmp_dir = tmpdir
        CONFIG.results.library_source = ResultsSource(
            file_1=Source(short_name="hsapiens", taxon_id=9606),
            file_2=Source(short_name="hsapiens", taxon_id=9606)
        )
        CONFIG.results.library_type = ResultsType(
            relationship=StatesTypeRelationship.split_mates,
        )
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.not_available,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_paired_isf(self, tmpdir):
        """Get read orientation for a paired-end stranded forward library."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_ISF_1
        CONFIG.args.path_2_processed = FILE_ORIENTATION_ISF_2
        CONFIG.results.library_source = ResultsSource(
            file_1=Source(short_name="hsapiens", taxon_id=9606),
            file_2=Source(short_name="hsapiens", taxon_id=9606),
        )
        CONFIG.results.library_type = ResultsType(
            relationship=StatesTypeRelationship.split_mates,
        )
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_forward,
            file_2=StatesOrientation.stranded_reverse,
            relationship=StatesOrientationRelationship.inward_stranded_forward,
        )

    def test_evaluate_paired_isr(self, tmpdir):
        """Get read orientation for a paired-end stranded reverse library."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_ISR_1
        CONFIG.args.path_2_processed = FILE_ORIENTATION_ISR_2
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_reverse,
            file_2=StatesOrientation.stranded_forward,
            relationship=StatesOrientationRelationship.inward_stranded_reverse,
        )

    def test_evaluate_paired_iu(self, tmpdir):
        """Get read orientation for a paired-end unstranded library."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_IU_1
        CONFIG.args.path_2_processed = FILE_ORIENTATION_IU_2
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.unstranded,
            file_2=StatesOrientation.unstranded,
            relationship=StatesOrientationRelationship.inward_unstranded,
        )

    def test_process_single_dummy_sam_file_file_problem(self, tmpdir):
        """Pass dummy aligned.out.sam file to simulate file
        problem."""
        CONFIG.args.tmp_dir = tmpdir
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        with pytest.raises(FileProblem):
            test_instance.process_single(FILE_EMPTY_ALIGNED_SAM)

    def test_process_paired_dummy_sam_file(self, tmpdir):
        """Pass dummy aligned.out.sam file to simulate file
        problem."""
        CONFIG.args.tmp_dir = tmpdir
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        with pytest.raises(FileProblem):
            test_instance.process_paired(FILE_EMPTY_ALIGNED_SAM)

    def test_process_paired_wrong_sam_file(self, tmpdir):
        """Pass bad_aligned.out.sam file to ensure correct
        paired file behaviour."""
        CONFIG.args.tmp_dir = tmpdir
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        assert test_instance.process_paired(FILE_BAD_ALIGNED_SAM) \
               == ResultsOrientation()

    def test_evaluate_paired_not_mates_unmapped(self, tmpdir):
        """Get read orientation for a paired-end library with no mappable
        reads.
        """
        CONFIG.args.path_1_processed = FILE_UNMAPPED_PAIRED_1
        CONFIG.args.path_2_processed = FILE_IDS_NOT_MATCH_2
        CONFIG.results.library_type = ResultsType(
                relationship=StatesTypeRelationship.not_available,
            )
        CONFIG.results.library_source = ResultsSource(
            file_1=Source(),
            file_2=Source(short_name="hsapiens", taxon_id=9606),
        )
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = GetOrientation(config=CONFIG,
                                       mapping=MAPPING)
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.not_available,
            file_2=StatesOrientation.stranded_reverse,
            relationship=StatesOrientationRelationship.not_available,
        )
