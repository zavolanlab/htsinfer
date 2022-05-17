"""Infer read orientation from sample data."""

import logging
from pathlib import Path
import tempfile
from typing import (Optional, Tuple)

from Bio import SeqIO  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
)
from htsinfer.models import (
    ResultsStats,
)

LOGGER = logging.getLogger(__name__)


class GetLibStats:
    """Determine library statitics of a single- or paired-end seguencing
    library.

    Args:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        tmp_dir: Path to directory where temporary output is written to.

    Attributes:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        tmp_dir: Path to directory where temporary output is written to.
    """
    def __init__(
        self,
        paths: Tuple[Path, Optional[Path]],
        tmp_dir: Path = Path(tempfile.gettempdir()) / 'tmp_htsinfer',
    ):
        """Class contructor."""
        self.paths = paths
        self.tmp_dir = tmp_dir

    def evaluate(self) -> ResultsStats:
        """Infer read statistics.

        Returns:
            Statistics results object.
        """
        stats = ResultsStats()

        # process file 1
        LOGGER.info(f"Obtaining statistics for file: {self.paths[0]}")
        stats.file_1.read_length.min, stats.file_1.read_length.max = (
            self.fastq_get_min_max_read_length(fastq=self.paths[0])
        )
        # process file 2
        LOGGER.info(f"Obtaining statistics for file: {self.paths[0]}")
        if self.paths[1] is not None:
            stats.file_2.read_length.min, stats.file_2.read_length.max = (
                self.fastq_get_min_max_read_length(fastq=self.paths[0])
            )

        return stats

    @staticmethod
    def fastq_get_min_max_read_length(fastq: Path) -> Tuple[int, int]:
        """Get number of records in a FASTQ file.

        Args:
            fastq: Path to FASTQ file.

        Returns:
            Tuple of minimum and maximum read lengths in input file.

        Raises:
            FileProblem: Could not process FASTQ file.
        """
        LOGGER.debug(
            f"Extracting minimum and maximum read length in: {fastq}")
        min_len: int = 1000000
        max_len: int = 0
        try:
            for record in SeqIO.parse(
                handle=fastq,
                format='fastq',
            ):
                min_len = min(min_len, len(record.seq))
                max_len = max(max_len, len(record.seq))
        except OSError as exc:
            raise FileProblem(
                f"Failed to process FASTQ file: {fastq}"
            ) from exc
        LOGGER.debug(
            f"Extracted minimum and maximum read length: {(min_len, max_len)}"
        )
        return (min_len, max_len)
