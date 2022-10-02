"""Unit tests for module ``get_read_orientation.py``."""

import pytest

from htsinfer.exceptions import FileProblem, StarProblem
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
    FILE_2000_RECORDS,
    FILE_DUMMY,
    FILE_EMPTY_ALIGNED_SAM,
    FILE_INVALID_TRANSCRIPTS,
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
    RaiseOSError,
    SubprocessError,
    SOURCE_HUMAN,
    SOURCE_FRUIT_FLY,
)


class TestGetOrientation:
    """Test ``GetOrientation`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        test_instance = GetOrientation(
            paths=(FILE_MATE_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
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
            library_source=ResultsSource(),
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
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmp_dir,
            threads_star=1,
            min_mapped_reads=20,
            min_fraction=0.75,

        )
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2
        assert test_instance.library_type == ResultsType()
        assert test_instance.library_source == ResultsSource()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS
        assert test_instance.tmp_dir == tmp_dir
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
            library_source=ResultsSource(),
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
            library_source=ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(),
            ),
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
            library_source=ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(),
            ),
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
            library_source=ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(),
            ),
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
            library_source=ResultsSource(),
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
            library_source=ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(short_name="hsapiens", taxon_id=9606),
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
            library_source=ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(short_name="hsapiens", taxon_id=9606),
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
            library_source=ResultsSource(
                file_1=Source(short_name="hsapiens", taxon_id=9606),
                file_2=Source(short_name="hsapiens", taxon_id=9606),
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

    def test_subset_transcripts_by_organism(self, tmpdir):
        """Get filtered orgainsm transcripts for different organisms."""
        library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY
            )
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, FILE_ORIENTATION_IU_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            library_source=library_source,
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.subset_transcripts_by_organism()
        filtered_organisms_transcripts = tmpdir / f"{library_source}.fasta"
        assert results == filtered_organisms_transcripts

    def test_subset_transcripts_by_organism_file_problem(self, tmpdir):
        """Pass dummy file as transcripts.fasta file to simulate
        file problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_DUMMY,
            tmp_dir=tmpdir
        )
        with pytest.raises(FileProblem):
            test_instance.subset_transcripts_by_organism()

    def test_subset_transcripts_by_organism_invalid_fasta(self, tmpdir):
        """Pass invalid transcripts.fasta file to simulate index error."""
        library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY
            )
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=library_source,
            transcripts_file=FILE_INVALID_TRANSCRIPTS,
            tmp_dir=tmpdir
        )
        results = test_instance.subset_transcripts_by_organism()
        filtered_organisms_transcripts = tmpdir / f"{library_source}.fasta"
        assert results == filtered_organisms_transcripts

    def test_get_fasta_size(self, tmpdir):
        """Get nucleotide statistics for filtererd transcripts
        with different organisms."""
        library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY,
            )
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, FILE_ORIENTATION_IU_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            library_source=library_source,
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        filtered_organisms_transcripts = \
            test_instance.subset_transcripts_by_organism()
        results = test_instance.get_fasta_size(filtered_organisms_transcripts)
        assert results == 249986

    def test_get_fasta_size_file_problem(self, tmpdir):
        """Pass dummy file as filtered_organisms_transcripts
        to simulate file problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.get_fasta_size(FILE_DUMMY)

    def test_get_star_index_string_size(self, tmpdir):
        """Get length of STAR SA pre-indexing string."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.get_star_index_string_size(249986)
        assert results == 7

    def test_evaluate_star_index_problem(self, monkeypatch, tmpdir):
        """Force raising exception to stimulate a star problem."""
        library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY,
            )
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, FILE_ORIENTATION_IU_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            ),
            library_source=library_source,
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        monkeypatch.setattr(
            'htsinfer.get_read_orientation.GetOrientation.create_star_index',
            lambda *args, **kwargs: StarProblem
        )
        with pytest.raises(StarProblem):
            test_instance.evaluate()

    def test_prepare_star_alignment_commands(self, tmpdir):
        """Get star alignment command."""
        test_instance = GetOrientation(
            paths=(FILE_2000_RECORDS, None),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.not_mates,
            ),
            library_source=ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=Source(),
            ),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        index_dir = tmpdir / 'index'
        file1_alignment_path = tmpdir / 'alignments/file_1'
        cmd = "STAR --alignIntronMax 1 --alignEndsType Local --runThreadN 1" \
            + " --genomeDir " + str(index_dir) + " --outFilterMultimapNmax " \
            + "50 --outSAMunmapped Within KeepPairs --readFilesIn " \
            + str(FILE_2000_RECORDS) + " --outFileNamePrefix " \
            + str(file1_alignment_path) + "/"
        results = test_instance.prepare_star_alignment_commands(
            index_dir=index_dir
            )
        assert ' '.join(list(results.values())[0]) == cmd

    def test_generate_star_alignments_problem(self, monkeypatch, tmpdir):
        """Force raising exception to simulate problem."""
        library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY,
            )
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, FILE_ORIENTATION_IU_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.not_mates,
            ),
            library_source=library_source,
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        sub_method_name = 'htsinfer.get_read_orientation.' + \
            'GetOrientation.generate_star_alignments'
        monkeypatch.setattr(
            sub_method_name,
            lambda *args, **kwargs: StarProblem,
        )
        with pytest.raises(FileProblem):
            test_instance.evaluate()

    def test_generate_star_alignments_dummy_cmd(self, tmpdir):
        """Pass dummy cmd to force simulate star problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        index_dir = tmpdir / 'index'
        file1_alignment_path = tmpdir / 'alignments/file_1'
        dummy_cmd = [
            'STAR', '--alignIntrnMax', '1',
            '--alignEndsType', 'Local', '--runThreadN', '1',
            "--genomeDir", f"{str(index_dir)}",
            ]
        cmds = {file1_alignment_path: dummy_cmd}
        with pytest.raises(StarProblem):
            test_instance.generate_star_alignments(cmds)

    def test_process_single_dummy_sam_file(self, tmpdir):
        """Pass dummy aligned.out.sam file to simulate file
        problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.process_single(FILE_EMPTY_ALIGNED_SAM)

    def test_process_paired_dummy_sam_file(self, tmpdir):
        """Pass dummy aligned.out.sam file to simulate file
        problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        with pytest.raises(FileProblem):
            test_instance.process_paired(FILE_EMPTY_ALIGNED_SAM)

    def test_create_star_index_star_problem(self, tmpdir):
        """Pass invalid transcripts path to simulate star problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        transcripts_path = tmpdir / 'invalid'
        with pytest.raises(StarProblem):
            test_instance.create_star_index(transcripts_path)

    def test_evaluate_paired_not_mates_unmapped(self, tmpdir):
        """Get read orientation for a paired-end library with no mappable
        reads.
        """
        test_instance = GetOrientation(
            paths=(FILE_UNMAPPED_PAIRED_1, FILE_UNMAPPED_PAIRED_2),
            library_type=ResultsType(
                relationship=StatesTypeRelationship.not_mates,
            ),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        results = test_instance.evaluate()
        assert results == ResultsOrientation(
            file_1=StatesOrientation.not_available,
            file_2=StatesOrientation.not_available,
            relationship=StatesOrientationRelationship.not_available,
        )

    def test_subset_transcripts_by_organism_cannot_write_file(
        self, monkeypatch, tmpdir
    ):
        """Force raising of ``OSError`` to simulate file problem."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_INVALID_TRANSCRIPTS,
            tmp_dir=tmpdir
        )
        monkeypatch.setattr(
            'Bio.SeqIO.write',
            RaiseOSError,
        )
        with pytest.raises(FileProblem):
            test_instance.subset_transcripts_by_organism()

    def test_generate_star_alignments_star_problem(self, monkeypatch, tmpdir):
        """Force raising of ``SubprocessError`` to simulate star probelm."""
        test_instance = GetOrientation(
            paths=(FILE_ORIENTATION_IU_1, None),
            library_type=ResultsType(),
            library_source=ResultsSource(),
            transcripts_file=FILE_TRANSCRIPTS,
            tmp_dir=tmpdir,
        )
        file1_alignment_path = tmpdir / 'alignments/file_1'
        dummy_cmd = [
            'STAR', '--alignIntrnMax', '1',
            '--alignEndsType', 'Local', '--runThreadN', '1',
            "--genomeDir",
            ]
        cmds = {file1_alignment_path: dummy_cmd}
        monkeypatch.setattr(
            'subprocess.run',
            lambda *args, **kwargs: SubprocessError(),
        )
        with pytest.raises(StarProblem):
            test_instance.generate_star_alignments(cmds)
