"""Infer read orientation from sample data."""

import logging
import os
import shutil
import subprocess as sp
import tempfile
from typing import Union
import zipfile as zp

from Bio import SeqIO  # type: ignore

LOGGER = logging.getLogger(__name__)


def infer(
    fasta: str,
    file_1: str,  # pylint: disable=unused-argument
    file_2: str = None,  # pylint: disable=unused-argument
    organism: Union[int, str] = "hsapiens",
    threads: int = 1,
) -> str:
    """Infers read orientation for single- or paired-ended sequencing libraries
    in FASTQ format.

    Note:
        File passed to `fasta` is expected to contain `|`-separated sequence
        identifier lines that contain an organism short name and a taxon
        identifier in the fourth and fifth columns, respectively. Example
        sequence identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`

    Args:
        fasta: File path to transcripts FASTA file; cf. note.
        file_1: File path to read/first mate library.
        file_2: File path to second mate library.
        organism: Source organism of the sequencing library; either an organism
            short name (string, e.g., `hsapiens`) or a taxon identifier
            (integer, e.g., `9606`).
        threads: Number of threads to run STAR.

    Returns:
        LIBTYPE string according to Salmon documentation, cf.
        https://salmon.readthedocs.io/en/latest/library_type.html
    """

    # create temporary directory in "../data" directory
    try:
        tmp_dir = tempfile.mkdtemp(dir=os.getcwd())
    except OSError as exc:
        raise OSError("Creation of temporary directory failed") from exc
    LOGGER.info(f"Created temporary directory '{tmp_dir}'")

    # Extracting fasta file
    with zp.ZipFile(fasta, "r") as zip_ref:
        zip_ref.extractall()

    # write FASTA file with transcripts of specified organism
    organism_transcripts = os.path.join(tmp_dir, f"{str(organism)}.fasta")
    subset_fasta_by_organism(
        fasta_in="transcripts.fasta",
        fasta_out=organism_transcripts,
        organism=organism,
    )
    LOGGER.info(
        f"Extracted sequence records for organism '{organism}' from FASTA "
        f"file '{fasta}' and wrote to FASTA file '{organism_transcripts}'"
    )

    # implement logic

    # Creating transcript index for mapping
    index_dir = star_index(
        tmp_dir=tmp_dir,
        organism_transcripts=organism_transcripts,
        organism=organism,
    )

    # delete temporary directory
    try:
        shutil.rmtree(tmp_dir)
    except OSError as exc:
        raise OSError("Deletion of temporary directory failed") from exc
    LOGGER.info(f"Deleted temporary directory '{tmp_dir}'")

    # return orientation string
    return "U"


def star_index(
    tmp_dir: str,
    organism_transcripts: str,
    organism: Union[int, str] = "hsapiens",
    threads: int = 1,
) -> str:
    """Prepares trascript index for mapping.
    Args:
        tmp_dir: Path to temprary directory.
        organism_transcripts: Extracted organism transcripts.
        organism: Source organism of the sequencing library; either an organism
            short name (string, e.g., `hsapiens`) or a taxon identifier
            (integer, e.g., `9606`).
        threads: Number of threads to run STAR.

    Returns:
        Path to the build index file.

    """
    index_dir = os.path.join(tmp_dir, "index")
    index_cmd = "STAR --runThreadN " + str(threads) + " --runMode " + \
        "genomeGenerate --genomeDir " + index_dir + " --genomeFastaFiles " + \
        organism_transcripts
    try:
        LOGGER.debug(f"Running '{organism}' transcript index")
        sp.run(index_cmd, shell=True, check=True)
        return index_dir
    except sp.CalledProcessError:
        LOGGER.error(
            "Error: running STAR."
        )
        return "E"


def subset_fasta_by_organism(
    fasta_in: str,
    fasta_out: str,
    organism: Union[int, str] = "hsapiens",
) -> None:
    """
    Writes FASTA records of the specified organism to file.

    Note:
        File passed to `fasta_in` is expected to contain `|`-separated sequence
        identifier lines that contain an organism short name and a taxon
        identifier in the fourth and fifth columns, respectively. Example
        sequence identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`

    Args:
        fasta_in: Path to FASTA input file.
        fasta_out: Path where FASTA output file shall be written to.
        organism: Either an organism short name (string, e.g., `hsapiens`) or a
            taxon identifier (integer, e.g., `9606`).

    Raises:
        OSError: Raised if `fasta_in` cannot be read from or `fasta_out`
            cannot be written to.
    """
    def yield_filtered_seqs():
        """Generator yielding sequence records for specified organism."""
        try:
            for record in SeqIO.parse(
                handle=fasta_in,
                format='fasta',
            ):
                try:
                    _, _, _, org_name, taxon_id = record.description.split("|")
                except ValueError:
                    continue
                if (
                    org_name == str(organism) or
                    taxon_id == str(organism)
                ):
                    yield record
        except OSError as exc:
            raise OSError(f"Could not open input file '{fasta_in}'") from exc

    try:
        SeqIO.write(
            sequences=yield_filtered_seqs(),
            handle=fasta_out,
            format='fasta',
        )
    except OSError as exc:
        raise OSError(f"Failed to write to FASTA file '{fasta_out}'") from exc
