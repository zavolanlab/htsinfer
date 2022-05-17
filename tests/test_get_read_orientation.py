"""Unit tests for module ``get_read_layout.py``."""

from htsinfer.get_read_orientation import GetOrientation
from htsinfer.models import (
    ResultsOrientation,
    ResultsType,
    StatesOrientation,
    StatesOrientationRelationship,
    StatesTypeRelationship,
)
from tests.utils import (
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
)


class TestGetOrientation:
    """Test ``GetOrientation`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        test_instance = GetOrientation(
            paths=(FILE_MATE_1, None),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
        )
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.library_type == ResultsType()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        test_instance = GetOrientation(
            paths=(FILE_MATE_1, FILE_MATE_2),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
        )
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2
        assert test_instance.library_type == ResultsType()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS

    def test_init_all(self, tmpdir):
        """Create instance with all available parameters."""
        tmp_dir = tmpdir
        test_instance = GetOrientation(
            paths=(FILE_MATE_1, FILE_MATE_2),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmp_dir,
            organism="hsapiens",
            threads_star=1,
            min_mapped_reads=20,
            min_fraction=0.75,

        )
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2
        assert test_instance.library_type == ResultsType()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS
        assert test_instance.tmp_dir == tmp_dir
        assert test_instance.organism == "hsapiens"
        assert test_instance.threads_star == 1
        assert test_instance.min_mapped_reads == 20
        assert test_instance.min_fraction == 0.75

    def test_evaluate_single_unmapped(self, tmpdir):
        """Get read orientation for a single-end library with no mappable
        reads.
        """
        test_instance = GetOrientation(
            paths=(FILE_UNMAPPED_SINGLE, None),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.not_available,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_single_sf(self, tmpdir):
        """Get read orientation for a single-end stranded forward library."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_SF, None),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_forward,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_single_sr(self, tmpdir):
        """Get read orientation for a single-end stranded reverse library."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_SR, None),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_reverse,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_single_u(self, tmpdir):
        """Get read orientation for a single-end unstranded library."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_U, None),
            library_type=ResultsType(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
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
        test_instance = GetOrientation(
            paths=(FILE_UNMAPPED_PAIRED_1, FILE_UNMAPPED_PAIRED_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.not_available,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_evaluate_paired_isf(self, tmpdir):
        """Get read orientation for a paired-end stranded forward library."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_ISF_1, FILE_ORIENTATION_ISF_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_forward,
            file_2=StatesOrientation.stranded_reverse,
            relationship=StatesOrientationRelationship.inward_stranded_forward,
        )

    def test_evaluate_paired_isr(self, tmpdir):
        """Get read orientation for a paired-end stranded reverse library."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_ISR_1, FILE_ORIENTATION_ISR_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.stranded_reverse,
            file_2=StatesOrientation.stranded_forward,
            relationship=StatesOrientationRelationship.inward_stranded_reverse,
        )

    def test_evaluate_paired_iu(self, tmpdir):
        """Get read orientation for a paired-end unstranded library."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, FILE_ORIENTATION_IU_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.unstranded,
            file_2=StatesOrientation.unstranded,
            relationship=StatesOrientationRelationship.inward_unstranded,
        )
