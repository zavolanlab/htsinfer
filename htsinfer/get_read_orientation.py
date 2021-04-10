"""Infer read orientation from sample data."""

from collections import defaultdict
import logging
from pathlib import Path
import subprocess as sp
import tempfile
from typing import (DefaultDict, Tuple, List, Optional)
import zipfile as zp

from Bio import SeqIO  # type: ignore
import pysam  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    MetadataWarning,
    StarProblem,
)
from htsinfer.models import ResultsOrientation

LOGGER = logging.getLogger(__name__)


class GetOrientation:
    """Determine the read orientation present in the FASTQ sequencing libraries.

    Note:
        File passed to `fasta` is expected to contain `|`-separated sequence
        identifier lines that contain an organism short name and a taxon
        identifier in the fourth and fifth columns, respectively. Example
        sequence identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`

    Args:
        fasta: File path to transcripts FASTA file.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        threads: Number of threads to run STAR.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        tmp_dir: Path to directory where temporary output is written to.

    Attrubutes:
        fasta: File path to transcripts FASTA file.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        threads: Number of threads to run STAR.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        tmp_dir: Path to directory where temporary output is written to.
        results: Results container for storing read orientation information for
            the provided files.
    """
    def __init__(
        self,
        fasta: Path,
        path_1: Path,
        path_2: Optional[Path] = None,
        threads: int = 1,
        organism: str = "hsapiens",
        tmp_dir: Path = Path(tempfile.gettempdir()),
    ):
        """Class contructor."""
        self.fasta: Path = fasta
        self.path_1: Path = path_1
        self.path_2: Optional[Path] = path_2
        self.threads = threads
        self.organism = organism
        self.tmp_dir = tmp_dir
        self.results: ResultsOrientation = ResultsOrientation()

    def evaluate(self) -> None:
        """Decides Orientation"""
        # Extracts transcripts fasta file
        try:
            with zp.ZipFile(self.fasta, "r") as zip_ref:
                zip_ref.extractall(self.tmp_dir)
        except (FileNotFoundError, zp.BadZipFile) as exc:
            self.results.file_1 = None
            self.results.file_2 = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        # process file 1
        LOGGER.debug(f"Processing file: '{self.path_1}'")
        orientation_file_1 = GetOutcomeType(
            fasta=self.fasta,
            path=self.path_1,
            threads=self.threads,
            organism=self.organism,
            tmp_dir=self.tmp_dir,
        )
        orientation_file_1.evaluate()
        self.results.file_1 = orientation_file_1.result
        LOGGER.debug(f"Orientation: {self.results.file_1}")

        # process file 2
        if self.path_2 is not None:
            LOGGER.debug(f"Processing file: '{self.path_2}'")
            orientation_file_2 = GetOutcomeType(
                fasta=self.fasta,
                path=self.path_2,
                threads=self.threads,
                organism=self.organism,
                tmp_dir=self.tmp_dir,
            )
            orientation_file_2.evaluate()
            self.results.file_2 = orientation_file_2.result
            LOGGER.debug(f"Orientation: {self.results.file_2}")


class GetOutcomeType():
    """Determine orientation for an individual FASTQ library.

    Args:
        fasta: File path to transcripts FASTA file.
        path: File path to read library.
        threads: Number of threads to run STAR.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        tmp_dir: Path to directory where temporary output is written to.

    Attributes:
        fasta: File path to transcripts FASTA file.
        path: File path to read library.
        threads: Number of threads to run STAR.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        tmp_dir: Path to directory where temporary output is written to.
        organism_transcripts: Extracted organism transcripts.
        alignment_count: Dictionary with alignment information for every read.
        result:  Orientation type for the FASTQ file.
    """
    def __init__(
        self,
        fasta: Path,
        path: Path,
        threads: int = 1,
        organism: str = "hsapiens",
        tmp_dir: Path = Path(tempfile.gettempdir()),
    ):
        """Class constructor."""
        self.fasta: Path = fasta
        self.path: Path = path
        self.threads = threads
        self.organism = organism
        self.tmp_dir: Path = tmp_dir
        self.organism_transcripts: Path = Path(tmp_dir) / \
            f"{str(organism)}.fasta"
        self.alignment_count: DefaultDict[str, List[int]] = defaultdict(
            lambda: [0, 0]
            )
        self.result: Optional[str] = None

    def evaluate(self) -> None:
        """Decides orientation type."""
        # Write FASTA file with transcripts of specified organism
        LOGGER.debug(
            f"Extracting sequence records for organism '{self.organism}' from "
            f"FASTA file '{self.fasta}' and writing to FASTA file "
            f"'{self.organism_transcripts}'"
        )
        try:
            _file = Path(self.tmp_dir) / Path(Path(self.fasta).name).stem
            self._subset_fasta_by_organism(
                fasta_in=_file,
                fasta_out=self.organism_transcripts,
            )
        except (OSError, ValueError) as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        try:
            self._star_index()
        except StarProblem as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        try:
            self._star_map()
        except StarProblem as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        # Extracts alignment from sam file
        self._extract_alignment()

        forward, reverse = self._count_alignment()
        self._decide_outcome(forward=forward, reverse=reverse)

    def _star_index(self) -> None:
        """Prepares trascript index for mapping."""
        LOGGER.debug(f"Running '{self.organism}' transcript index")
        _dir = Path(self.tmp_dir) / "index"
        index_cmd = "STAR --runThreadN " + str(self.threads) + " --runMode" + \
            " genomeGenerate --genomeDir " + str(_dir) + \
            " --genomeFastaFiles " + str(self.organism_transcripts) + \
            " --genomeSAindexNbases 7"
        result = sp.run(index_cmd, shell=True, capture_output=True, text=True)
        LOGGER.debug(result.stderr)
        if result.returncode == 0:
            pass
        else:
            raise StarProblem("Failed to run STAR index")

    def _star_map(self) -> None:
        """Maps reads to transcripts."""
        LOGGER.debug(f"Mapping reads to transcripts for '{self.fasta}'")
        _dir = Path(self.tmp_dir) / "alignment/"  # CHECK
        map_cmd = "STAR --runThreadN " + str(self.threads) + \
            " --genomeDir " + str(self.tmp_dir) + "/index --readFilesIn " + \
            str(self.path) + " --alignIntronMax 1 --alignEndsType Local " + \
            "--outFileNamePrefix " + str(_dir) + "/"
        result = sp.run(map_cmd, shell=True, capture_output=True, text=True)
        LOGGER.debug(result.stderr)
        if result.returncode == 0:
            pass
        else:
            raise StarProblem("Failed to run STAR mapping job")

    def _extract_alignment(self) -> None:
        """Extracts alignment configuration from Aligned.out.sam file."""
        _file = Path(self.tmp_dir) / "alignment" / "Aligned.out.sam"
        LOGGER.debug("Extracting alignment configuration")
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
            self.result = None
            raise MetadataWarning(
                f"Failed to extract alignment configuration for '{self.path}'"
                ) from exc

    def _count_alignment(self) -> Tuple[int, int]:
        """Decides the outcomes of all alignments of a given read into one
        outcome.

        Returns:
            A tuple of two items representing the count of alignemnts for sense
            and anti-sense.
        """
        forward_f: float = 0
        reverse_f: float = 0
        total_f: int = 0

        if not self.alignment_count:
            return (0, 0)
        for _fwd, _rev in self.alignment_count.values():
            forward_f += _fwd / (_fwd + _rev)
            reverse_f += _rev / (_fwd + _rev)
            total_f += 1

        return (
                round((forward_f/total_f)*100), round((reverse_f/total_f)*100)
            )

    def _decide_outcome(
        self,
        forward: int,
        reverse: int,
    ) -> None:
        """Uses the read alignment configuration counts to decide the orientation type.
        Args:
            forward: Count alignment percentage towards forward.
            reverse: Count alignment percentage towards reverse.

        Note:
            stranded, reads in forward configuration : SF
            stranded, reads in reverse-complemented configuration : SR
            unstranded : U
        """
        if forward == 0 and reverse == 0:
            self.result = None
        elif forward >= 80:
            self.result = "SF"
        elif forward <= 20:
            self.result = "SR"
        elif 40 <= forward <= 60:
            self.result = "U"
        else:
            self.result = None

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
