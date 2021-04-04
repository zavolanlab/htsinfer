"""Infer read orientation from sample data."""

from collections import defaultdict
import logging
import os
import shutil
import subprocess as sp
import tempfile
from typing import (DefaultDict, Tuple, List, Union)
import zipfile as zp

from Bio import SeqIO  # type: ignore
import pysam  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
)
from htsinfer.models import ResultsOrientation

LOGGER = logging.getLogger(__name__)


class GetOrientation:
    """Determine the read orientation present in the FASTQ sequencing libraries.

    Args:
        fasta: File path to transcripts FASTA file.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        threads: Number of threads to run STAR.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
    
    Attrubutes:
        fasta: File path to transcripts FASTA file.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        threads: Number of threads to run STAR.
        organism: Source organism of the sequencing library, if provided:
            will not not be inferred by the application.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        results: Results container for storing read orientation information for
            the provided files.
    """


