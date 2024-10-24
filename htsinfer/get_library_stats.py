"""Infer read orientation from sample data."""

import logging
from pathlib import Path
from collections import Counter
import statistics
from typing import Tuple

from Bio import SeqIO  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
)
from htsinfer.models import (
    ResultsStats,
    Config,
)

LOGGER = logging.getLogger(__name__)


class GetLibStats:
    """
    Determine library statistics of a single- or paired-end sequencing library.

    Args:
        config: Container class for all arguments used in inference
                and results produced by the class.

    Attributes:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        tmp_dir: Path to directory where temporary output is written to.
    """

    def __init__(
        self,
        config: Config,
    ):
        """Class contructor."""
        self.paths = (config.args.path_1_processed,
                      config.args.path_2_processed)
        self.tmp_dir = config.args.tmp_dir / 'tmp_htsinfer'

    def evaluate(self) -> ResultsStats:
        """Infer read statistics.

        Returns:
            Statistics results object.
        """
        stats = ResultsStats()

        # process file 1
        LOGGER.info(f"Obtaining statistics for file: {self.paths[0]}")
        (stats.file_1.read_length.min,
         stats.file_1.read_length.max,
         stats.file_1.read_length.mean,
         stats.file_1.read_length.median,
         stats.file_1.read_length.mode) = (
            self.fastq_get_stats_read_length(fastq=self.paths[0])
        )
        # process file 2
        LOGGER.info(f"Obtaining statistics for file: {self.paths[1]}")
        if self.paths[1] is not None:
            (stats.file_2.read_length.min,
             stats.file_2.read_length.max,
             stats.file_2.read_length.mean,
             stats.file_2.read_length.median,
             stats.file_2.read_length.mode) = (
                self.fastq_get_stats_read_length(fastq=self.paths[1])
            )

        return stats

    @staticmethod
    def fastq_get_stats_read_length(fastq: Path) -> Tuple[
        int, int, float, int, int
            ]:
        """Get number of records in a FASTQ file.

        Args:
            fastq: Path to FASTQ file.

        Returns:
            Tuple of minimum and maximum read lengths in input file.

        Raises:
            FileProblem: Could not process FASTQ file.
        """
        LOGGER.debug(
            f"Extracting read length statistics in: {fastq}"
        )
        min_len: int = 1000000
        max_len: int = 0
        total_lengths: int = 0
        lengths: list = []
        length_counter: Counter = Counter()

        try:
            for record in SeqIO.parse(
                handle=fastq,
                format='fastq',
            ):
                seq_length = len(record.seq)
                min_len = min(min_len, seq_length)
                max_len = max(max_len, seq_length)
                total_lengths += seq_length
                lengths.append(seq_length)
                length_counter[seq_length] += 1
        except OSError as exc:
            raise FileProblem(
                f"Failed to process FASTQ file: {fastq}"
            ) from exc

        mean_len = round(total_lengths / len(lengths), 2)
        median_len = int(statistics.median(lengths))
        mode_len = length_counter.most_common(1)[0][0]

        LOGGER.debug(
            f"Extracted read length statistics: "
            f"{(min_len, max_len, mean_len, median_len, mode_len)}"
        )
        return min_len, max_len, mean_len, median_len, mode_len
