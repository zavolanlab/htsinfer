"""Infer read orientation from sample data."""

from typing import Union


def infer(
    file_1: str,  # pylint: disable=unused-argument
    file_2: str = None,  # pylint: disable=unused-argument
    organism: Union[int, str] = "hsapiens",  # pylint: disable=unused-argument
) -> str:
    """Infers read orientation for single- or paired-ended sequencing libraries
    in FASTQ format.

    Args:
        file_1: File path to read/first mate library.
        file_2: File path to second mate library.
        organism: Source organism of the sequencing library; either a short
            name (string, e.g., `hsapiens`) or a taxon identifier (integer,
            e.g., `9606`).

    Returns:
        LIBTYPE string according to Salmon documentation, cf.
        https://salmon.readthedocs.io/en/latest/library_type.html
    """
    # implement logic
    return "U"
