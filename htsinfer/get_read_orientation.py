"""Infer read orientation from sample data."""

from collections import defaultdict
import logging
import math
from pathlib import Path
import subprocess as sp
from typing import (Any, DefaultDict, Dict, List)

from Bio import SeqIO  # type: ignore
import pysam  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    StarProblem,
)
from htsinfer.models import (
    ResultsOrientation,
    StatesOrientation,
    StatesOrientationRelationship,
    StatesTypeRelationship,
    Config,
)

LOGGER = logging.getLogger(__name__)


class GetOrientation:
    """Determine library strandedness and relative read orientation of a
    single- or paired-end seguencing library.

    Args:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        library_type: ResultsType object with library type and mate
            relationship.
        library_source: ResultsSource object with source information on each
            library file.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format.
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        source: Source (organism, tissue, etc.) of the sequencing library.
        min_mapped_reads: Minimum number of mapped reads for deeming the
            read orientation result reliable.
        min_fraction: Minimum fraction of mapped reads required to be
            consistent with a given read orientation state in order for that
            orientation to be reported. Must be above 0.5.
        mate_relationship: Type/mate relationship between the provided files.

    Attributes:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        library_type: ResultsType object with library type and mate
            relationship.
        library_source: ResultsSource object with source information on each
            library file.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format.
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        source: Source (organism, tissue, etc.) of the sequencing library.
        min_mapped_reads: Minimum number of mapped reads for deeming the
            read orientation result reliable.
        min_fraction: Minimum fraction of mapped reads required to be
            consistent with a given read orientation state in order for that
            orientation to be reported. Must be above 0.5.
        mate_relationship: Type/mate relationship between the provided files.
    """
    def __init__(
        self,
        config: Config,
    ):
        """Class contructor."""
        self.paths = (config.args.path_1_processed,
                      config.args.path_2_processed)
        self.library_type = config.results.library_type
        self.library_source = config.results.library_source
        self.transcripts_file = config.args.t_file_processed
        self.tmp_dir = config.args.tmp_dir
        self.threads_star = config.args.threads
        self.min_mapped_reads = config.args.read_orientation_min_mapped_reads
        self.min_fraction = config.args.read_orientation_min_fraction

    def evaluate(self) -> ResultsOrientation:
        """Infer read orientation.

        Returns:
            Orientation results object.
        """
        orientation = ResultsOrientation()

        # get transcripts for current organims
        transcripts = self.subset_transcripts_by_organism()
        ref_size = self.get_fasta_size(fasta=transcripts)
        index_string_size = self.get_star_index_string_size(ref_size=ref_size)

        # generate STAR alignments
        index_dir = self.create_star_index(
            fasta=transcripts,
            index_string_size=index_string_size,
        )
        star_cmds = self.prepare_star_alignment_commands(index_dir=index_dir)
        self.generate_star_alignments(commands=star_cmds)

        # process alignments
        star_dirs = [e for e in star_cmds if e is not None]
        orientation = self.process_alignments(star_dirs=star_dirs)

        return orientation

    def subset_transcripts_by_organism(self) -> Path:
        """Filter FASTA file of transcripts by current sources.

        The filtered file contains records from the indicated sources.
            Typically, this is one source. However, for if two input files
            were supplied that are originating from different sources (i.e.,
            not from a valid paired-ended library), it may be from two
            different sources. If no source is supplied (because it could
            not be inferred), no filtering is done.

        Returns:
            Path to filtered FASTA file.

        Raises:
            FileProblem: Could not open input/output FASTA file for
                reading/writing.
        """
        LOGGER.debug(f"Subsetting transcripts for: {self.library_source}")

        outfile = self.tmp_dir / f"{self.library_source}.fasta"

        def yield_filtered_seqs():
            """Generator yielding sequence records for specified sources.

            Yields:
                Next FASTA sequence record of the specified sources.

            Raises: Could not process input FASTA file.
            """
            sources = []
            if self.library_source.file_1.short_name is not None:
                sources.append(self.library_source.file_1.short_name)
            if self.library_source.file_2.short_name is not None:
                sources.append(self.library_source.file_2.short_name)
            try:
                for record in SeqIO.parse(
                    handle=self.transcripts_file,
                    format='fasta',
                ):
                    try:
                        org_name = record.description.split("|")[3]
                    except ValueError:
                        continue
                    if org_name in sources or len(sources) == 0:
                        yield record

            except OSError as exc:
                raise FileProblem(
                    f"Could not process file '{self.transcripts_file}'"
                ) from exc

        try:
            SeqIO.write(
                sequences=yield_filtered_seqs(),
                handle=outfile,
                format='fasta',
            )
        except OSError as exc:
            raise FileProblem(
                f"Failed to write to FASTA file '{outfile}'"
            ) from exc

        LOGGER.debug(f"Filtered transcripts file: {outfile}")
        return outfile

    @staticmethod
    def get_fasta_size(fasta: Path) -> int:
        """Get size of FASTA file in total nucleotides.

        Args:
            fasta: Path to FASTA file.

        Returns:
            Total number of nucleotides of all records.

        Raises:
            FileProblem: Could not open FASTA file for reading.
        """
        nucleotides: int = 0

        try:
            for record in SeqIO.parse(
                handle=fasta,
                format='fasta',
            ):
                nucleotides += len(record.seq)

        except OSError as exc:
            raise FileProblem(
                f"Could not process file: {fasta}"
            ) from exc

        LOGGER.debug(f"Size of reference: {nucleotides}")
        return nucleotides

    @staticmethod
    def get_star_index_string_size(ref_size: int) -> int:
        """Get length of STAR SA pre-indexing string.

        Cf.
        https://github.com/alexdobin/STAR/blob/51b64d4fafb7586459b8a61303e40beceeead8c0/doc/STARmanual.pdf

        Args:
            ref_size: Size of genome/transcriptome reference in nucleotides.

        Returns:
            Size (in nucleotides) of SA pre-indexing string.
        """
        index_string_size = min(
            14,
            int(math.floor(math.log2(ref_size) / 2 - 1))
        )
        LOGGER.debug(f"STAR SA pre-indexing string size: {index_string_size}")
        return index_string_size

    def create_star_index(
        self,
        fasta: Path,
        index_string_size: int = 5,
    ) -> Path:
        """Prepare STAR index.

        Args:
            fasta: Path to FASTA file of sequence records to create index from.
            index_string_size: Size of SA pre-indexing string, in nucleotides.

        Returns:
            Path to directory containing STAR index.

        Raises:
            StarProblem: STAR index could not be created.
        """
        LOGGER.debug(f"Creating STAR index for: {fasta}")

        index_dir: Path = Path(self.tmp_dir) / "index"

        cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeSAindexNbases", f"{str(index_string_size)}",
            "--runThreadN", f"{str(self.threads_star)}",
            "--genomeDir", f"{str(index_dir)}",
            "--genomeFastaFiles", f"{str(fasta)}",
        ]

        result = sp.run(
            cmd,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            LOGGER.error(result.stderr)
            raise StarProblem("Failed to create STAR index")

        LOGGER.debug(f"STAR index created: {index_dir}")
        return index_dir

    def prepare_star_alignment_commands(
        self,
        index_dir: Path,
    ) -> Dict[Path, List[str]]:
        """Prepare STAR alignment commands.

        Args:
            index_dir: Path to directory containing STAR index.

        Returns:
            Dictionary of output paths and corresponding STAR commands.
        """
        LOGGER.debug("Preparing STAR commands...")

        # helper function for compiling individual command
        def build_star_command(
            read_files: List[str],
            out_dir: str,
        ) -> List[str]:
            """Compile an individual STAR alignment command.

            Args:
                read_files: List of read file paths.
                out_dir: STAR output directory.

            Returns:
                STAR command list.
            """
            cmd_base: List[str] = [
                "STAR",
                "--alignIntronMax", "1",
                "--alignEndsType", "Local",
                "--runThreadN", f"{str(self.threads_star)}",
                "--genomeDir", f"{str(index_dir)}",
                "--outFilterMultimapNmax", "50",
                "--outSAMunmapped", "Within", "KeepPairs",
            ]
            cmd: List[str] = cmd_base[:]
            cmd.append("--readFilesIn")
            cmd.extend(read_files)
            cmd.append("--outFileNamePrefix")
            cmd.append(out_dir)
            return cmd

        out_dir_base: Path = Path(self.tmp_dir) / "alignments"
        commands: Dict = {}

        # create command for pairend-ended libraries
        if (
            self.library_type.relationship
            == StatesTypeRelationship.split_mates
        ):
            out_dir = out_dir_base / "paired"
            commands[out_dir] = build_star_command(
                read_files=[str(path) for path in self.paths],
                out_dir=f"{str(out_dir)}/",
            )
        # create commands for single-ended libraries
        else:
            out_dir = out_dir_base / "file_1"
            commands[out_dir] = build_star_command(
                read_files=[str(self.paths[0])],
                out_dir=f"{str(out_dir)}/",
            )
            # run two commands in case there is a second file provided that is
            # not a mate of the first one
            out_dir = out_dir_base / "file_2"
            if self.paths[1] is not None:
                commands[out_dir] = build_star_command(
                    read_files=[str(self.paths[1])],
                    out_dir=f"{str(out_dir)}/",
                )

        return commands

    @staticmethod
    def generate_star_alignments(commands: Dict[Path, List[str]]) -> None:
        """Align reads to index with STAR.

        Args:
            commands: Dictionary of output paths and corresponding STAR
                commands.

        Raises:
            StarProblem: Generating alignments failed.
        """
        LOGGER.debug("Aligning reads with STAR...")

        # execute commands
        for out_dir, cmd in commands.items():
            try:
                result = sp.run(
                    cmd,
                    capture_output=True,
                    text=True,
                )
                if result.returncode != 0:
                    LOGGER.error(result.stderr)
                    raise StarProblem(
                        "Failed to generate STAR alignments for command: "
                        f"{cmd}"
                    )
            except sp.CalledProcessError as exc:
                raise StarProblem(
                    f"Failed to generate STAR alignments for command: {cmd}"
                ) from exc
            LOGGER.debug(f"Written STAR output to directory: {out_dir}")

    def process_alignments(
        self,
        star_dirs: List[Path],
    ) -> ResultsOrientation:
        """Determine read orientation of one or two single-ended or one
        paired-end sequencing library.

        Args:
            star_dirs: List of one or two paths to STAR output directories.

        Returns:
            Read orientation state of library or libraries.
        """
        results: ResultsOrientation = ResultsOrientation()
        paths = [path / 'Aligned.out.sam' for path in star_dirs]
        if len(star_dirs) == 1:
            if star_dirs[0].name == "paired":
                results = self.process_paired(sam=paths[0])
            if star_dirs[0].name == "file_1":
                results.file_1 = self.process_single(sam=paths[0])
        elif len(star_dirs) == 2:
            results.file_1 = self.process_single(sam=paths[0])
            results.file_2 = self.process_single(sam=paths[1])
        return results

    def process_single(
        self,
        sam: Path,
    ) -> StatesOrientation:
        """Determine read orientation of a single-ended sequencing library.

        Args:
            sam: Path to SAM file.

        Returns:
            Read orientation state of library.
        """
        LOGGER.debug(f"Processing SAM file: '{sam}'")

        # parsing alignments
        states: DefaultDict[str, List[StatesOrientation]] = defaultdict(
            lambda: []
        )
        try:
            with pysam.AlignmentFile(str(sam), "r") as _file:
                for record in _file.fetch():

                    # ensure that "unmapped" flag is not set and query name
                    # is set
                    if (
                        not record.flag & (1 << 2) and
                        isinstance(record.query_name, str)
                    ):

                        # check if reverse complemented
                        if record.flag & (1 << 4):
                            states[record.query_name].append(
                                StatesOrientation.stranded_reverse
                            )
                        elif not record.flag & (1 << 4):
                            states[record.query_name].append(
                                StatesOrientation.stranded_forward
                            )

        except OSError as exc:
            raise FileProblem(
                f"Failed to open SAM file: '{sam}'"
            ) from exc

        LOGGER.debug("Deciding read orientation...")
        reads = len(states)
        fractions = [
            self.get_frequencies(*state) for state in states.values()
        ]
        fractions_all_states = {
            key: val / reads for (key, val) in
            self.sum_dicts(*fractions).items()
        }
        orientation: StatesOrientation = StatesOrientation.not_available
        if len(fractions_all_states) == 0:
            return orientation
        most_common_state = max(
            fractions_all_states,
            key=fractions_all_states.get  # type: ignore
        )
        fraction_most_common_state = max(fractions_all_states.values())
        if reads >= self.min_mapped_reads:
            if fraction_most_common_state > self.min_fraction:
                orientation = most_common_state
            else:
                orientation = StatesOrientation.unstranded

        # write log messages and return result
        LOGGER.debug(
            f"Required number of mapped reads pairs: {self.min_mapped_reads}"
        )
        LOGGER.debug(f"Number of reads mapped: {reads}")
        LOGGER.debug(f"Fraction of states: {fractions_all_states}")
        LOGGER.debug(f"Orientation: {orientation}")
        return orientation

    def process_paired(  # pylint: disable=R0912,R0915
        self,
        sam: Path,
    ) -> ResultsOrientation:
        """Determine read orientation of a paired-ended sequencing library.

        Args:
            sam: Path to SAM file.

        Returns:
            Read orientation state of each mate and orientation state
                relationship of library.
        """
        LOGGER.debug(f"Processing SAM file: '{sam}'")

        # parsing alignments
        states: DefaultDict[
            str, List[StatesOrientationRelationship]
        ] = defaultdict(lambda: [])
        try:
            with pysam.AlignmentFile(str(sam), "r") as _file:
                for record_1 in _file.fetch():

                    # ensure both mates are mapped and part of mate pair
                    if not (
                        record_1.flag & (1 << 0) and
                        record_1.flag & (1 << 1) and
                        not record_1.flag & (1 << 2) and
                        not record_1.flag & (1 << 3)
                    ):
                        continue

                    # check which alignment is first mate
                    record_2 = next(_file)
                    if (
                        record_1.flag & (1 << 6) and
                        record_2.flag & (1 << 7)
                    ):
                        mate_1 = record_1
                        mate_2 = record_2
                    elif (
                        record_2.flag & (1 << 6) and
                        record_1.flag & (1 << 7)
                    ):
                        mate_1 = record_2
                        mate_2 = record_1
                    else:
                        continue

                    # ensure that query name is present
                    if not isinstance(mate_1.query_name, str):
                        continue

                    # check orientation: forward / inward
                    if (
                        not mate_1.flag & (1 << 4) and
                        mate_1.pos < mate_2.pos  # type: ignore
                    ):
                        states[mate_1.query_name].append(
                            StatesOrientationRelationship.
                            inward_stranded_forward
                        )

                    # check orientation: reverse / inward
                    elif (
                        mate_1.flag & (1 << 4) and
                        mate_1.pos > mate_2.pos  # type: ignore
                    ):
                        states[mate_1.query_name].append(
                            StatesOrientationRelationship.
                            inward_stranded_reverse
                        )

        except StopIteration:
            pass

        except OSError as exc:
            raise FileProblem(
                f"Failed to open SAM file: '{sam}'"
            ) from exc

        # deciding read orientation
        LOGGER.debug("Deciding read orientation...")
        reads = len(states)
        fractions = [
            self.get_frequencies(*state) for state in states.values()
        ]
        fractions_all_states = {
            key: val / reads for (key, val) in
            self.sum_dicts(*fractions).items()
        }
        orientation: ResultsOrientation = ResultsOrientation()
        if len(fractions_all_states) == 0:
            return orientation
        most_common_state = max(
            fractions_all_states,
            key=fractions_all_states.get  # type: ignore
        )
        fraction_most_common_state = max(fractions_all_states.values())
        if reads >= self.min_mapped_reads:
            if fraction_most_common_state > self.min_fraction:
                orientation.relationship = most_common_state
                assert orientation.relationship is not None
                if orientation.relationship.value[-2:] == "SF":
                    orientation.file_1 = StatesOrientation("SF")
                    orientation.file_2 = StatesOrientation("SR")
                else:
                    orientation.file_1 = StatesOrientation("SR")
                    orientation.file_2 = StatesOrientation("SF")
            else:
                orientation.relationship = (
                    StatesOrientationRelationship.inward_unstranded
                )
                orientation.file_1 = StatesOrientation.unstranded
                orientation.file_2 = StatesOrientation.unstranded

        # write log messages and return result
        LOGGER.debug(
            f"Required number of mapped read pairs: {self.min_mapped_reads}"
        )
        LOGGER.debug(f"Number of reads mapped: {reads}")
        LOGGER.debug(f"Fraction of states: {fractions_all_states}")
        LOGGER.debug(f"Orientation: {orientation}")
        return orientation

    @staticmethod
    def get_frequencies(*items: Any) -> Dict[Any, float]:
        """Get frequencies of arguments as fractions of the number of all
        arguments.

        Args:
            *items: Items to get frequencies for.

        Returns:
            Dictionary of arguments and their frequencies.
        """
        counts: DefaultDict[Any, int] = defaultdict(int)
        length: int = len(items)
        for item in items:
            counts[item] += 1
        fractions: Dict[Any, float] = {
            k: v/length for (k, v) in counts.items()
        }
        return fractions

    @staticmethod
    def sum_dicts(*dicts: Dict[Any, float]) -> Dict[Any, float]:
        """Sum of dictionaries with numeric values.

        Args:
            *dicts: Dictionaries to sum up.

        Returns:
            Dictionary with union of keys of input dictionaries and all values
            added up.
        """
        result: DefaultDict[Any, float] = defaultdict(int)
        for dct in dicts:
            for key, num in dct.items():
                result[key] += num
        return dict(result)
