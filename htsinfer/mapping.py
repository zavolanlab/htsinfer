"""Mapping FASTQ's and managing the outputs of STAR."""

import logging
import math
from pathlib import Path
import subprocess as sp
from typing import (Dict, List)

from Bio import SeqIO  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    StarProblem,
)
from htsinfer.models import (
    Config,
    StatesTypeRelationship,
)

LOGGER = logging.getLogger(__name__)


class Mapping:
    """Map FASTQ file/s and manage outputs.

    Args:
        path: Path to FASTQ file.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.

    Raise:
        FileProblem: The input file could not be parsed or the output file
            could not be written.
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
        self.mapped = False
        self.star_dirs: List[Path] = []

    def evaluate(self):
        """Align FASTQ files to reference transcripts with STAR."""

        # get transcripts for current organims
        transcripts = self.subset_transcripts_by_organism()
        ref_size = self.get_fasta_size(fasta=transcripts)
        index_string_size = self.get_star_index_string_size(ref_size=ref_size)
        chr_bin_bits = self.get_star_chr_bin_bits(ref_size=ref_size,
                                                  transcripts=transcripts)

        # generate STAR alignments
        index_dir = self.create_star_index(
            fasta=transcripts,
            index_string_size=index_string_size,
            chr_bin_bits=chr_bin_bits,
        )
        star_cmds = self.prepare_star_alignment_commands(index_dir=index_dir)
        self.generate_star_alignments(commands=star_cmds)
        # process alignments
        self.star_dirs = [e for e in star_cmds if e is not None]
        self.mapped = True

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

        outfile = self.tmp_dir / "transcripts_subset.fasta"

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
                    except (ValueError, IndexError):
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

    @staticmethod
    def get_star_chr_bin_bits(ref_size: int, transcripts: Path) -> int:
        """Get size of bins for STAR genome storage.

        Args:
            ref_size: Size of genome/transcriptome reference in nucleotides.
            transcripts: Path to filtered FASTA transcripts file.

        Returns:
            Number of bins for genome storage.
        """
        n_ref: int = 0

        for _ in SeqIO.parse(
                handle=transcripts,
                format='fasta',
        ):
            n_ref += 1

        chr_bin_bits = min(
            18,
            int(round(math.log2(ref_size / n_ref)))
        )
        LOGGER.debug("STAR size of bins for genome storage: %s", chr_bin_bits)
        return chr_bin_bits

    def create_star_index(
            self,
            fasta: Path,
            chr_bin_bits: int = 18,
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

        # solves the macOS issue with STAR
        index_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeSAindexNbases", f"{str(index_string_size)}",
            "--genomeChrBinNbits", f"{str(chr_bin_bits)}",
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

        Input FASTQ files are assumed to be sorted according to reference names
        or coordinates, the order of input reads is kept with the option
        "PairedKeepInputOrder", no additional sorting of aligned reads is done.

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
                "--outSAMorder", "PairedKeepInputOrder",
                "--outSAMunmapped", "Within", "KeepPairs",
            ]
            cmd: List[str] = cmd_base[:]
            cmd.append("--readFilesIn")
            cmd.extend(read_files)
            cmd.append("--outFileNamePrefix")
            cmd.append(out_dir)

            # solves the macOS issue with STAR
            Path(out_dir).mkdir(parents=True, exist_ok=True)

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
                    check=True,
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
