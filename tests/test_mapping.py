"""Unit tests for module ``mapping.py``."""

import pytest

from htsinfer.exceptions import (
    FileProblem,
    StarProblem,
)
from htsinfer.mapping import Mapping
from htsinfer.models import (
    ResultsSource,
    ResultsType,
    Source,
    StatesTypeRelationship,
)
from tests.utils import (
    FILE_2000_RECORDS,
    FILE_DUMMY,
    FILE_INVALID_TRANSCRIPTS,
    FILE_MATE_1,
    FILE_MATE_2,
    FILE_ORIENTATION_IU_1,
    FILE_ORIENTATION_IU_2,
    FILE_TRANSCRIPTS,
    CONFIG,
    MAPPING,
    RaiseError,
    SubprocessError,
    SOURCE_HUMAN,
    SOURCE_FRUIT_FLY,
)


class TestMapping:
    """Test ``Mapping`` class."""

    def test_init_required(self):
        """Create instance with required parameters."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = None
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.results.library_type = ResultsType()
        test_instance = Mapping(config=CONFIG)
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.library_type == ResultsType()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS

    def test_init_required_paired(self):
        """Create instance with required parameters for paired-end library."""
        CONFIG.args.path_1_processed = FILE_MATE_1
        CONFIG.args.path_2_processed = FILE_MATE_2
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        test_instance = Mapping(config=CONFIG)
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
        test_instance = Mapping(config=CONFIG)
        assert test_instance.paths[0] == FILE_MATE_1
        assert test_instance.paths[1] == FILE_MATE_2
        assert test_instance.library_type == ResultsType()
        assert test_instance.library_source == ResultsSource()
        assert test_instance.transcripts_file == FILE_TRANSCRIPTS
        assert test_instance.tmp_dir == tmpdir

    def test_subset_transcripts_by_organism(self, tmpdir):
        """Get filtered orgainsm transcripts for different organisms."""
        CONFIG.results.library_type = ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            )
        CONFIG.results.library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY
            )
        CONFIG.args.tmp_dir = tmpdir
        MAPPING.mapped = False
        test_instance = Mapping(config=CONFIG)
        results = test_instance.subset_transcripts_by_organism()
        filtered_organisms_transcripts = \
            tmpdir / "transcripts_subset.fasta"
        assert results == filtered_organisms_transcripts

    def test_subset_transcripts_by_organism_file_problem(self, tmpdir):
        """Pass dummy file as transcripts.fasta file to simulate
        file problem."""
        CONFIG.args.path_2_processed = None
        CONFIG.results.library_type = ResultsType()
        CONFIG.results.library_source = ResultsSource()
        CONFIG.args.t_file_processed = FILE_DUMMY
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)

        with pytest.raises(FileProblem):
            test_instance.subset_transcripts_by_organism()

    def test_subset_transcripts_by_organism_invalid_fasta(self, tmpdir):
        """Pass invalid transcripts.fasta file to simulate index error."""
        CONFIG.results.library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY
        )
        CONFIG.args.t_file_processed = FILE_INVALID_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        results = test_instance.subset_transcripts_by_organism()
        filtered_organisms_transcripts = \
            tmpdir / "transcripts_subset.fasta"
        assert results == filtered_organisms_transcripts

    def test_get_fasta_size(self, tmpdir):
        """Get nucleotide statistics for filtererd transcripts
        with different organisms."""
        CONFIG.results.library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY,
            )
        CONFIG.args.path_2_processed = FILE_ORIENTATION_IU_2
        CONFIG.results.library_type = ResultsType(
            relationship=StatesTypeRelationship.split_mates,
        )
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        filtered_organisms_transcripts = \
            test_instance.subset_transcripts_by_organism()
        results = test_instance.get_fasta_size(filtered_organisms_transcripts)
        assert results == 249986

    def test_get_fasta_size_file_problem(self, tmpdir):
        """Pass dummy file as filtered_organisms_transcripts
        to simulate file problem."""
        CONFIG.args.path_2_processed = None
        CONFIG.results.library_type = ResultsType()
        CONFIG.results.library_source = ResultsSource()
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        with pytest.raises(FileProblem):
            test_instance.get_fasta_size(FILE_DUMMY)

    def test_get_star_index_string_size(self, tmpdir):
        """Get length of STAR SA pre-indexing string."""
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        results = test_instance.get_star_index_string_size(249986)
        assert results == 7

    def test_evaluate_star_index_problem(self, monkeypatch, tmpdir):
        """Force raising exception to stimulate a star problem."""
        CONFIG.results.library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=SOURCE_FRUIT_FLY,
            )
        CONFIG.args.path_2_processed = FILE_ORIENTATION_IU_2
        CONFIG.results.library_type = ResultsType(
                relationship=StatesTypeRelationship.split_mates,
            )
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        monkeypatch.setattr(
            'htsinfer.mapping.Mapping.create_star_index',
            lambda *args, **kwargs: StarProblem
        )
        with pytest.raises(StarProblem):
            test_instance.evaluate()

    def test_prepare_star_alignment_commands(self, tmpdir):
        """Get star alignment command."""
        CONFIG.args.path_1_processed = FILE_2000_RECORDS
        CONFIG.args.path_2_processed = None
        CONFIG.results.library_type = ResultsType(
            relationship=StatesTypeRelationship.not_mates,
        )
        CONFIG.results.library_source = ResultsSource(
                file_1=SOURCE_HUMAN,
                file_2=Source(),
            )
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        index_dir = tmpdir / 'index'
        file1_alignment_path = tmpdir / 'alignments/file_1'
        cmd = "STAR --alignIntronMax 1 --alignEndsType Local --runThreadN 1" \
            + " --genomeDir " + str(index_dir) + " --outFilterMultimapNmax " \
            + "50 --outSAMorder PairedKeepInputOrder " \
            + "--outSAMunmapped Within KeepPairs --readFilesIn " \
            + str(FILE_2000_RECORDS) + " --outFileNamePrefix " \
            + str(file1_alignment_path) + "/"
        results = test_instance.prepare_star_alignment_commands(
            index_dir=index_dir
            )
        assert ' '.join(list(results.values())[0]) == cmd

    def test_generate_star_alignments_dummy_cmd(self, tmpdir):
        """Pass dummy cmd to force simulate star problem."""
        CONFIG.args.path_2_processed = None
        CONFIG.results.library_type = ResultsType()
        CONFIG.results.library_source = ResultsSource()
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
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

    def test_create_star_index_star_problem(self, tmpdir):
        """Pass invalid transcripts path to simulate star problem."""
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        transcripts_path = tmpdir / 'invalid'
        with pytest.raises(StarProblem):
            test_instance.create_star_index(transcripts_path)

    def test_subset_transcripts_by_organism_cannot_write_file(
        self, monkeypatch, tmpdir
    ):
        """Force raising of ``OSError`` to simulate file problem."""
        CONFIG.args.path_1_processed = FILE_ORIENTATION_IU_1
        CONFIG.args.path_2_processed = None
        CONFIG.results.library_source = ResultsSource()
        CONFIG.results.library_type = ResultsType()
        CONFIG.args.t_file_processed = FILE_INVALID_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
        monkeypatch.setattr(
            'Bio.SeqIO.write',
            RaiseError(exc=OSError),
        )
        with pytest.raises(FileProblem):
            test_instance.subset_transcripts_by_organism()

    def test_generate_star_alignments_star_problem(self, monkeypatch, tmpdir):
        """Force raising of ``SubprocessError`` to simulate star probelm."""
        CONFIG.args.t_file_processed = FILE_TRANSCRIPTS
        CONFIG.args.tmp_dir = tmpdir
        test_instance = Mapping(config=CONFIG)
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
