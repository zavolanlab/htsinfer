"""Infer read layout from sample data."""

from typing import Union


def infer(
    file_1: str,  # pylint: disable=unused-argument
    file_2: str = None,  # pylint: disable=unused-argument
    organism: Union[int, str] = "hsapiens",  # pylint: disable=unused-argument
) -> str:
    """Infers read layout for single- or paired-ended sequencing libraries in
    FASTQ format.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.

    Returns:
        LIBTYPE string according to Salmon documentation, cf.
        https://salmon.readthedocs.io/en/latest/library_type.html
    """
    # implement logic
    return "U"
