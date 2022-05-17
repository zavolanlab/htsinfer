"""FASTQ subsetting, extraction and validation."""

from functools import partial
import gzip
import logging
from pathlib import Path

from Bio import SeqIO  # type: ignore

from htsinfer.exceptions import FileProblem

LOGGER = logging.getLogger(__name__)


class SubsetFastq:
    """Subset, uncompress and validate a FASTQ file.

    Args:
        path: Path to FASTQ file.
        out_dir: Path to directory where output is written to.
        records: Number of input file records to process; set to `0` to
            process all records.

    Attributes:
        path: Path to FASTQ file.
        out_dir: Path to directory where output is written to.
        records: Number of input file records to process.
        out_path: Path for uncompressed, filtered `path` file.
        n_processed: Total number of processed records.

    Raise:
        FileProblem: The input file could not be parsed or the output file
            could not be written.
    """
    def __init__(
        self,
        path: Path,
        out_dir: Path = Path.cwd() / 'results_htsinfer',
        records: int = 0,
    ):
        """Class constructor."""
        self.path = path
        self.out_dir = out_dir
        self.records = records
        self.out_path: Path = self.path
        self.n_processed = 0

    def process(self):
        """Uncompress, subset and validate files."""
        LOGGER.debug(f"Processing file: {self.path}")
        if self.path.suffix == ".gz":
            name = self.path.stem
            _open = partial(gzip.open)
        else:
            name = self.path.name
            _open = open
        self.out_path = self.out_dir / name
        try:
            with open(self.out_path, 'wt', encoding="utf-8") as out_file:
                with _open(self.path, 'rt') as in_file:
                    seq_iter = SeqIO.parse(
                        handle=in_file,
                        format='fastq',
                    )
                    for record in seq_iter:
                        self.n_processed += SeqIO.write(
                            sequences=record,
                            handle=out_file,
                            format='fastq',
                        )
                        if (
                            self.records and
                            self.n_processed >= self.records
                        ):
                            break
                if not self.n_processed:
                    raise FileProblem(f"File is empty: {self.path}")
        except Exception as exc:
            raise FileProblem(exc) from exc
        LOGGER.debug(
            f"Written {self.n_processed} records to: {self.out_path}"
        )
