"""Infer read orientation from sample data."""

import logging
import os
import shutil
from tempfile import mkdtemp
from typing import Union

LOGGER = logging.getLogger(__name__)


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

    # create temporary directory in "../data" directory
    try:
        tmp_dir = mkdtemp(dir=os.getcwd())
    except OSError as exc:
        raise OSError("Creation of temporary directory failed") from exc
    LOGGER.info(f"Created temporary directory '{tmp_dir}'")

    # implement logic

    # delete temporary directory
    try:
        shutil.rmtree(tmp_dir)
    except OSError as exc:
        raise OSError("Deletion of temporary directory failed") from exc
    LOGGER.info(f"Deleted temporary directory '{tmp_dir}'")

    # return orientation string
    return "U"
