"""Infer read orientation from sample data."""

from collections import defaultdict
import logging
from pathlib import Path
import subprocess as sp
import tempfile
from typing import (DefaultDict, Tuple, List, Optional)

from Bio import SeqIO  # type: ignore
import pysam  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    StarProblem,
)
from htsinfer.models import (
    ResultsOrientation,
    ResultsType,
    StatesOrientation,
    StatesOrientationRelationship,
    StatesTypeRelationship,
)

LOGGER = logging.getLogger(__name__)


class GetOrientation:
    """Determine library strandedness and relative read orientation of a
    single- or paired-end seguencing library.

    Args:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        library_type: Library type and mate relationship.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        organism: Source organism of the sequencing library.
        min_mapped_reads: Minimum number of mapped reads for deeming the
            read orientation result reliable.
        fraction_range: Size of the range of the fraction of mapped reads that
            are consistent with one of the outcomes 'stranded-forward',
            'stranded-reverse' and 'unstranded'. Must be at least zero and at
            most one third.
        mate_relationship: Type/mate relationship between the provided files.

    Attributes:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        library_type: Library type and mate relationship.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        organism: Source organism of the sequencing library.
        min_mapped_reads: Minimum number of mapped reads for deeming the
            read orientation result reliable.
        fraction_range: Size of the range of the fraction of mapped reads that
            are consistent with one of the outcomes 'stranded-forward',
            'stranded-reverse' and 'unstranded'. Must be at least zero and at
            most one third.
        mate_relationship: Type/mate relationship between the provided files.
    """
    def __init__(
        self,
        paths: Tuple[Path, Optional[Path]],
        library_type: ResultsType,
        transcripts_file: Path = (
            Path(__file__).parent.parent.absolute() /
            "data/transcripts.fasta.gz"
        ),
        tmp_dir: Path = Path(tempfile.gettempdir()) / 'tmp_htsinfer',
        threads_star: int = 1,
        organism: str = "hsapiens",
        min_mapped_reads: int = 20,
        fraction_range: float = 0.2,
    ):
        """Class contructor."""
        self.paths = paths
        self.library_type = library_type
        self.transcripts_file = transcripts_file
        self.tmp_dir = tmp_dir
        self.threads_star = threads_star
        self.organism = organism
        self.min_mapped_reads = min_mapped_reads
        self.fraction_range = fraction_range

    def evaluate(self) -> ResultsOrientation:
        """Infer read orientation.

        Returns:
            Orientation results object.
        """
        orientation = ResultsOrientation()

        # get transcripts for current organims
        transcripts: Path = self.subset_transcripts_by_organism()

        # create STAR index
        index_dir: Path = self.create_star_index(fasta=transcripts)

        # process first library
        orientation.file_1 = self.process_library(
            path=self.paths[0],
            index_dir=index_dir,
        )

        # process second library, if available
        if self.paths[1] is not None:
            orientation.file_2 = self.process_library(
                path=self.paths[1],
                index_dir=index_dir,
            )
            orientation.relationship = self.decide_relationship(
                orientation_1=orientation.file_1,
                orientation_2=orientation.file_2,
            )

        return orientation

    def subset_transcripts_by_organism(self) -> Path:
        """Filter FASTA file of transcripts by current organism.

        Returns:
            Path to filtered FASTA file.

        Raises:
            FileProblem: Could not open input/output FASTA file for
                reading/writing.
        """
        LOGGER.debug(f"Subsetting transcripts for '{self.organism}'...")

        outfile = self.tmp_dir / f"{self.organism}.fasta"

        def yield_filtered_seqs():
            """Generator yielding sequence records for specified organism.

            Yields:
                Next FASTA sequence record of the specified organism.
            """
            try:
                for record in SeqIO.parse(
                    handle=self.transcripts_file,
                    format='fasta',
                ):
                    try:
                        org_name = record.description.split("|")[3]
                    except ValueError:
                        continue
                    if org_name == self.organism:
                        yield record

            except OSError as exc:
                raise FileProblem(
                    f"Could not open input file '{self.transcripts_file}'"
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

    def create_star_index(
        self,
        fasta: Path,
    ) -> Path:
        """Prepare STAR index.

        Args:
            fasta: Path to FASTA file of sequence records to create index from.

        Returns:
            Path to directory containing STAR index.

        Raises:
            StarProblem: STAR index could not be created.
        """
        LOGGER.debug(f"Creating STAR index for '{fasta}'...")

        index_dir: Path = Path(self.tmp_dir) / "index"

        cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeSAindexNbases", "7",
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

    def generate_star_alignments(
        self,
        path: Path,
        index_dir: Path,
    ) -> Path:
        """Align reads to index with STAR.

        Args:
            path: File path to read library.
            index_dir: Path to directory containing STAR index.

        Returns:
            STAR output directory.
        """
        LOGGER.debug(f"Aligning reads for '{path}'...")

        out_dir: Path = Path(self.tmp_dir) / "alignment"

        cmd = [
            "STAR",
            "--alignIntronMax", "1",
            "--alignEndsType", "EndToEnd",
            "--runThreadN", f"{str(self.threads_star)}",
            "--genomeDir", f"{str(index_dir)}",
            "--readFilesIn", f"{str(path)}",
            "--outFileNamePrefix", f"{str(out_dir)}/",
        ]

        result = sp.run(
            cmd,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            LOGGER.error(result.stderr)
            raise StarProblem("Failed to generate STAR alignments")

        LOGGER.debug(f"STAR alignemnts generated: {out_dir}")
        return out_dir

    @staticmethod
    def count_alignment_configurations(
        sam: Path,
    ) -> Tuple[int, float]:
        """Get alignment configuration from SAM file.

        Returns:
            Tuple of
            (1) Total number of "mapped" reads;
            (2) Total number of reads that were not reverse complemented before
                being aligned to the reference. A floating point number,
                because for reads with alignments to multiple locations, the
                individual configurations for each alignment are averaged. For
                example, if a read has four equivalent alignments, and in one
                of these the read first needed to be reverse complemented
                before aligning, the overall contribution of this individual
                read to the total number of reads is 0.75.

        Raises:
            FileProblem: SAM file could not be opened.
        """
        LOGGER.debug(f"Extracting alignment configurations for '{sam}'...")

        # create dictionary of alignment configurations for each read
        # keys: read identifiers
        # values: list of configurations for each alignment, each
        #   represented by either 0 or 1, with,
        #   0: alignment to sense strand
        #   1: alignment to antisense strand
        configs: DefaultDict[str, List[int]] = defaultdict(lambda: [])
        try:
            with pysam.AlignmentFile(sam, "r") as _file:
                for record in _file.fetch():
                    # ensure that "unmapped" flag is not set
                    if not record.flag & (1 << 2):
                        # check that "reverse complemented" flag is not set
                        if not record.flag & (1 << 4):
                            configs[record.query_name].append(1)
                        else:
                            configs[record.query_name].append(0)
        except OSError as exc:
            raise FileProblem(
                f"Failed to open SAM file '{sam}'"
            ) from exc
        count_reads = len(configs)

        # calculate total number of reads consistent with aligning to the
        # 'sense'/'+' strand
        sum_fractions: float = 0
        for config in configs.values():
            sum_fractions += sum(config) / len(config)

        LOGGER.debug(f"Number of mapped reads: {count_reads}")
        LOGGER.debug(
            f"Number of reads aligning to 'sense'/'+' strand: {sum_fractions}"
        )
        return (len(configs), sum_fractions)

    def decide_orientation(
        self,
        count_reads: int,
        count_sense: float,
    ) -> StatesOrientation:
        """Decide orientation of library based on alignments to sense strand.

        Args:
            count_reads: Total number of "mapped" reads.
            count_sense: Total number of reads that were not reverse
                complemented before being aligned to the reference. A floating
                point number, because for reads with alignments to multiple
                locations, the individual configurations for each alignment are
                averaged. For example, if a read has four equivalent
                alignments, and in one of these the read first needed to be
                reverse complemented before aligning, the overall contribution
                of this individual read to the total number of reads is 0.75.
        """
        LOGGER.debug("Deciding orientation...")

        orientation: StatesOrientation = StatesOrientation.not_available

        try:
            fraction_sense = count_sense / count_reads
        # ensure we have at least one mapped read
        except ZeroDivisionError:
            LOGGER.warning(
                "Could not decide read orientation of file. No reads mapped."
            )
            return orientation

        # ensure enough records were considered
        if not count_reads > self.min_mapped_reads:
            LOGGER.warning(
                "Not enough reads evaluated to decide read orientation of "
                "file. Try increasing the number of records to consider or "
                "lower the minimum number of mapped reads."
            )

        # check for forward orientation
        if fraction_sense >= 1.0 - self.fraction_range:
            orientation = StatesOrientation.stranded_forward

        # check for reverse orientation
        elif fraction_sense <= 0.0 + self.fraction_range:
            orientation = StatesOrientation.stranded_reverse

        # check for unstranded
        elif (
            0.5 - (self.fraction_range / 2) <=
            fraction_sense <=
            0.5 + (self.fraction_range / 2)
        ):
            orientation = StatesOrientation.unstranded

        LOGGER.debug(f"Read orientation: {orientation}")
        return orientation

    def process_library(
        self,
        path: Path,
        index_dir: Path,
    ) -> StatesOrientation:
        """Determine read orientation of a single seguencing library file.

        Args:
            path: File path to read library.
            index_dir: Path to directory containing STAR index.

        Returns:
            Read orientation state of library.
        """
        alignments_dir = self.generate_star_alignments(
            path=path,
            index_dir=index_dir,
        )

        count_reads, count_sense = self.count_alignment_configurations(
            sam=alignments_dir / "Aligned.out.sam",
        )

        return self.decide_orientation(
            count_reads=count_reads,
            count_sense=count_sense,
        )

    def decide_relationship(
        self,
        orientation_1: StatesOrientation,
        orientation_2: StatesOrientation,
    ) -> StatesOrientationRelationship:
        """Determine relative configuration of mates.

        Args:
            orientation_1: Read orientation of first mate.
            orientation_2: Read orientation of second mate.

        Returns:
            Relative orientation of mate pairs.
        """
        relationship: StatesOrientationRelationship = (
            StatesOrientationRelationship.not_available
        )

        if not (
            self.library_type.relationship
            == StatesTypeRelationship.split_mates
        ):
            LOGGER.warning(
                "Could not decide orientation relationship because the "
                "sequence files are not mates of a paired-end library"
            )

        # unavailable
        elif (
            orientation_1 == StatesOrientation.not_available or
            orientation_2 == StatesOrientation.not_available
        ):
            LOGGER.warning(
                "Could not decide orientation relationship because "
                "orientation of mates could not be decided"
            )

        # inward_stranded_forward
        elif (
            orientation_1 == StatesOrientation.stranded_forward and
            orientation_2 == StatesOrientation.stranded_reverse
        ):
            relationship = (
                StatesOrientationRelationship.inward_stranded_forward
            )

        # inward_stranded_reverse
        elif (
            orientation_1 == StatesOrientation.stranded_reverse and
            orientation_2 == StatesOrientation.stranded_forward
        ):
            relationship = (
                StatesOrientationRelationship.inward_stranded_reverse
            )

        # TO IMPLEMENT
        # inward_unstranded
        # matching_stranded_forward
        # matching_stranded_reverse
        # matching_unstranded
        # outward_stranded_forward
        # outward_stranded_reverse

        return relationship
