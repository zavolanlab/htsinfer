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
    MetadataWarning,
    StarProblem,
)
from htsinfer.models import (
    ResultsOrientation,
    StatesOrientation,
    StatesOrientationRelationship,
)

LOGGER = logging.getLogger(__name__)


class GetOrientation:
    """Determine library strandedness and relative read orientation of a
    single- or paired-end seguencing library.

    Args:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        organism: Source organism of the sequencing library.
        range: Size of the fraction ranges that determine one of three
            outcomes: stranded-forward (1 - `range`), stranded-reverse (0 +
            `range`) and unstranded (0.5 +- `range` / 2). Value must be at
            least 0.05 and at most one third.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        organism: Source organism of the sequencing library.
        range: Size of the fraction ranges that determine one of three
            outcomes: stranded-forward (1 - `range`), stranded-reverse (0 +
            `range`) and unstranded (0.5 +- `range` / 2). Value must be at
            least 0.05 and at most one third.
        results: Results container for storing read orientation information.
    """
    def __init__(
        self,
        path_1: Path,
        path_2: Optional[Path] = None,
        transcripts_file: Path = (
            Path(__file__).parent.parent.absolute() /
            "data/transcripts.fasta.gz"
        ),
        tmp_dir: Path = Path(tempfile.gettempdir()) / 'tmp_htsinfer',
        threads_star: int = 1,
        organism: str = "hsapiens",
        range: float = 0.2,
    ):
        """Class contructor."""
        self.path_1: Path = Path(path_1)
        self.path_2: Optional[Path] = None if path_2 is None else Path(path_2)
        self.transcripts: Path = Path(transcripts_file)
        self.tmp_dir = Path(tmp_dir)
        self.threads_star = threads_star
        self.organism = organism
        self.range = range
        self.results: ResultsOrientation = ResultsOrientation()

    def evaluate(self) -> None:
        """Decides Orientation"""

        # process file 1
        LOGGER.debug(f"Processing file: '{self.path_1}'")
        orientation_file_1 = GetOrientationType(
            path=self.path_1,
            transcripts_file=self.transcripts,
            tmp_dir=self.tmp_dir,
            threads_star=self.threads_star,
            organism=self.organism,
            range=self.range,
        )
        orientation_file_1.evaluate()
        self.results.file_1 = orientation_file_1.result
        LOGGER.debug(f"Orientation: {self.results.file_1}")

        # process file 2
        if self.path_2 is not None:
            LOGGER.debug(f"Processing file: '{self.path_2}'")
            orientation_file_2 = GetOrientationType(
                path=self.path_2,
                transcripts_file=self.transcripts,
                tmp_dir=self.tmp_dir,
                threads_star=self.threads_star,
                organism=self.organism,
                range=self.range,
            )
            orientation_file_2.evaluate()
            self.results.file_2 = orientation_file_2.result
            LOGGER.debug(f"Orientation: {self.results.file_2}")

            # determine relationship
            # TODO: implement
            LOGGER.debug(f"Orientation paired-end: {self.results}")
            self.results.relationship = (
                StatesOrientationRelationship.not_available
            )


class GetOrientationType():
    """Determine library strandedness and read orientation of a single
    seguencing library file.

    Args:
        path: File path to read library.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        range: Size of the fraction ranges that determine one of three
            outcomes: stranded-forward (1 - `range`), stranded-reverse (0 +
            `range`) and unstranded (0.5 +- `range` / 2). Value must be at
            least 0.05 and at most one third.

    Attributes:
        path: File path to read library.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        organism: Source organism of the sequencing library.
        range: Size of the fraction ranges that determine one of three
            outcomes: stranded-forward (1 - `range`), stranded-reverse (0 +
            `range`) and unstranded (0.5 +- `range` / 2). Value must be at
            least 0.05 and at most one third.
        organism_transcripts: Extracted organism transcripts.
        alignment_count: Dictionary with alignment information for every read.
        result: Orientation type for the FASTQ file.
    """
    def __init__(
        self,
        path: Path,
        transcripts_file: Path,
        tmp_dir: Path = Path(tempfile.gettempdir()),
        threads_star: int = 1,
        organism: str = "hsapiens",
        range: float = 0.2,
    ):
        """Class constructor."""
        self.path: Path = path
        self.transcripts_file: Path = transcripts_file
        self.tmp_dir: Path = tmp_dir
        self.threads_star = threads_star
        self.organism = organism
        self.range = range
        self.organism_transcripts: Path = Path(tmp_dir) / \
            f"{str(organism)}.fasta"
        self.alignment_count: DefaultDict[str, List[int]] = defaultdict(
            lambda: [0, 0]
            )
        self.result: StatesOrientation = StatesOrientation.not_available

    def evaluate(self) -> None:
        """Decides orientation type."""
        # write FASTA file with transcripts of specified organism
        LOGGER.debug(
            f"Extracting sequence records for organism '{self.organism}' from "
            f"FASTA file '{self.transcripts_file}' and writing to FASTA file "
            f"'{self.organism_transcripts}'"
        )
        try:
            _file = Path(self.tmp_dir) / Path(
                Path(self.transcripts_file).name
            ).stem
            self._subset_fasta_by_organism(
                fasta_in=_file,
                fasta_out=self.organism_transcripts,
            )
        except (OSError, ValueError) as exc:
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        try:
            self._star_index()
        except StarProblem as exc:
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        try:
            self._star_map()
        except StarProblem as exc:
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        # extract and count alignments
        self._extract_alignments()
        fraction_forward, _ = self._count_alignments()

        # decide outcome
        self._decide_outcome(
            fraction_forward=fraction_forward,
            range=self.range,
        )

    def _subset_fasta_by_organism(
        self,
        fasta_in: Path,
        fasta_out: Path,
    ) -> None:
        """Writes FASTA records of the specified organism to file.

        Args:
            fasta_in: Path to FASTA input file.
            fasta_out: Path where FASTA output file shall be written to.
            organism: Source organism of the sequencing library, if provided:
                will not not be inferred by the application.
        """
        def _yield_filtered_seqs():
            """Generator yielding sequence records for specified organism."""
            try:
                for record in SeqIO.parse(
                    handle=fasta_in,
                    format='fasta',
                ):
                    try:
                        _, _, _, org_name, _ = record.description.split("|")
                    except ValueError:
                        continue
                    if org_name == self.organism:
                        yield record
            except OSError as exc:
                raise OSError(
                    f"Could not open input file '{fasta_in}'"
                ) from exc

        try:
            SeqIO.write(
                sequences=_yield_filtered_seqs(),
                handle=fasta_out,
                format='fasta',
            )
        except OSError as exc:
            raise OSError(
                f"Failed to write to FASTA file '{fasta_out}'"
            ) from exc

    def _star_index(self) -> None:
        """Prepares trascript index for mapping."""
        LOGGER.debug(f"Creating STAR index for '{self.transcripts_file}'")
        _dir = Path(self.tmp_dir) / "index"
        index_cmd = (
            "STAR "
            "--runMode=genomeGenerate "
            "--genomeSAindexNbases 7 "
            f"--runThreadN={str(self.threads_star)} "
            f"--genomeDir={str(_dir)} "
            f"--genomeFastaFiles={str(self.organism_transcripts)} "
        )
        result = sp.run(
            index_cmd,
            shell=True,
            capture_output=True,
            text=True,
        )
        LOGGER.debug(result.stderr)
        if result.returncode != 0:
            raise StarProblem("Failed to create STAR index")

    def _star_map(self) -> None:
        """Maps reads to transcripts."""
        LOGGER.debug(f"Mapping reads file '{self.path}' to transcript index")
        out_dir = str(Path(self.tmp_dir) / "alignment/")
        map_cmd = (
            "STAR "
            "--alignIntronMax=1 "
            "--alignEndsType=Local "
            f"--runThreadN={str(self.threads_star)} "
            f"--genomeDir={str(self.tmp_dir)}/index "
            f"--readFilesIn={str(self.path)} "
            f"--outFileNamePrefix={str(out_dir)}/ "
        )
        result = sp.run(
            map_cmd,
            shell=True,
            capture_output=True,
            text=True,
        )
        LOGGER.debug(result.stderr)
        if result.returncode != 0:
            raise StarProblem("Failed to run STAR mapping job")

    def _extract_alignments(self) -> None:
        """Extracts alignment configuration from Aligned.out.sam file."""
        _file = Path(self.tmp_dir) / "alignment" / "Aligned.out.sam"
        LOGGER.debug("Extracting alignment configurations")
        try:
            sam_file = pysam.AlignmentFile(_file, "r")
            for line in sam_file.fetch():
                flag = line.flag
                read = line.query_name
                unmapped: bool = flag & (1 << 2)
                if not unmapped:
                    if flag & (1 << 4):
                        self.alignment_count[read][1] += 1
                    else:
                        self.alignment_count[read][0] += 1
            sam_file.close()
        except OSError as exc:
            self.result = StatesOrientation.not_available
            raise MetadataWarning(
                f"Failed to extract alignment configuration for '{self.path}'"
            ) from exc

    def _count_alignments(self) -> Tuple[float, float]:
        """Calculate fractions of read alignments compatible with a forward
        and reverse read orientation.

        Returns:
            A tuple of two items representing the fractions of alignments
            compatible with a forward and reverse read orientation,
            respectively.
        """
        LOGGER.debug("Counting alignment configurations")
        forward: float = 0
        reverse: float = 0
        total: int = 0
        if not self.alignment_count:
            return (0, 0)
        for _fwd, _rev in self.alignment_count.values():
            forward += _fwd / (_fwd + _rev)
            reverse += _rev / (_fwd + _rev)
            total += 1
        return forward/total, reverse/total

    def _decide_outcome(
        self,
        fraction_forward: float,
        range: float = 0.2,
    ) -> None:
        """Decide orientation type based on read alignment configuration.

        Args:
            fraction_forward: Fraction of read alignments compatible with
                forward orientation.
            range: Size of the fraction ranges that determine one of three
                outcomes: stranded-forward (1 - `range`), stranded-reverse
                (0 + `range`) and unstranded (0.5 +- `range` / 2). Value must
                be at least 0.05 and at most one third.
        """
        LOGGER.debug("Deciding outcome")
        if not 0.05 <= range <= 1/3:
            raise ValueError(
                f"Parameter range ({range}) outside the permitted boundaries: "
                "must be between 0.05 and 1/3"
            )
        if fraction_forward >= 1.0 - range:
            self.result = StatesOrientation.stranded_forward
        elif fraction_forward <= 0.0 + range:
            self.result = StatesOrientation.stranded_reverse
        elif 0.5 - (range/2) <= fraction_forward <= 0.5 + (range/2):
            self.result = StatesOrientation.unstranded
