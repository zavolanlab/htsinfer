"""Infer library source from sample data."""

from collections import defaultdict
import logging
from pathlib import Path
import subprocess as sp
import tempfile
from typing import (DefaultDict, Dict, Optional, Tuple)

import pandas as pd  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    KallistoProblem,
)
from htsinfer.models import (
    ResultsSource,
    Source,
)
from htsinfer.utils import (
    convert_dict_to_df,
    validate_top_score,
)


LOGGER = logging.getLogger(__name__)


class GetLibSource:
    """Determine the source of FASTQ sequencing of a single- or paired-end
    seguencing library.

    Args:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        min_match_pct: Minimum percentage of reads that are consistent with a
            given source in order for it to be considered as the to be
            considered the library's source.
        min_freq_ratio: Minimum frequency ratio between the first and second
            most frequent source in order for the former to be considered the
            library's source.

    Attrubutes:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format. Expected to contain `|`-separated sequence identifier
            lines that contain an organism short name and a taxon identifier in
            the fourth and fifth columns, respectively. Example sequence
            identifier: `rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        min_match_pct: Minimum percentage of reads that are consistent with a
            given source in order for it to be considered as the to be
            considered the library's source.
        min_freq_ratio: Minimum frequency ratio between the first and second
            most frequent source in order for the former to be considered the
            library's source.
    """
    def __init__(
        self,
        paths: Tuple[Path, Optional[Path]],
        transcripts_file: Path,
        out_dir: Path = Path(
            __file__
        ).parents[2].absolute() / 'results_htsinfer',
        tmp_dir: Path = Path(tempfile.gettempdir()) / 'tmp_htsinfer',
        min_match_pct: float = 2,
        min_freq_ratio: float = 2,
    ):
        """Class contructor."""
        self.paths = paths
        self.transcripts_file = transcripts_file
        self.out_dir = out_dir
        self.tmp_dir = tmp_dir
        self.min_match_pct = min_match_pct
        self.min_freq_ratio = min_freq_ratio

    def evaluate(self) -> ResultsSource:
        """Infer read source.

        Returns:
            Source results object.
        """
        source = ResultsSource()
        index = self.create_kallisto_index()
        source.file_1 = self.get_source(
            fastq=self.paths[0],
            index=index,
        )
        if self.paths[1] is not None:
            source.file_2 = self.get_source(
                fastq=self.paths[0],
                index=index,
            )
        return source

    def create_kallisto_index(self) -> Path:
        """Build Kallisto index from FASTA file of target sequences.

        Returns:
            Path to Kallisto index.

        Raises:
            KallistoProblem: Kallisto index could not be created.
        """
        LOGGER.debug(f"Creating kallisto index for: {self.transcripts_file}")

        index = self.tmp_dir / "kallisto.idx"

        cmd = [
            "kallisto",
            "index",
            "--index", f"{str(index)}",
            "--make-unique",
            f"{str(self.transcripts_file)}",
        ]

        result = sp.run(
            cmd,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            LOGGER.error(result.stderr)
            raise KallistoProblem("Failed to run Kallisto index")

        LOGGER.debug(f"Kallisto index created: {index}")
        return index

    def get_source(
        self,
        fastq: Path,
        index: Path,
    ) -> Source:
        """Determine source of a single sequencing library file.

        Args:
            fastq: Path to FASTQ file.
            index: Path to Kallisto index.

        Returns:
            Source of library file.
        """
        source: Source = Source()

        # run quantification
        kallisto_dir = self.run_kallisto_quantification(
            fastq=fastq,
            index=index,
        )

        # extract counts
        tpm_pct = self.get_source_expression(
            kallisto_dir=kallisto_dir,
        )

        # convert counts to data frame
        sources_df = convert_dict_to_df(
            dic=tpm_pct,
            col_headers=('Source', 'Count Percentage'),
            sort=True,
            sort_by=1,
            sort_ascending=False,
        )

        # write data frame (in JSON) to file
        name = (
            Path(self.out_dir) / f"library_source_{fastq.name}.json"
        )
        LOGGER.debug(f"Writing results to file: {name}")
        sources_df.to_json(
            name,
            orient='split',
            index=False,
            indent=True,
        )

        # validate results
        if validate_top_score(
            vector=sources_df['Count Percentage'].to_list(),
            min_value=self.min_match_pct,
            min_ratio=self.min_freq_ratio,
            rev_sorted=True,
            accept_zero=True,
        ):
            LOGGER.warning(sources_df.iloc[0]['Source'])
            source.short_name, source.taxon_id = sources_df.iloc[0]['Source']

        return source

    def run_kallisto_quantification(
        self,
        fastq: Path,
        index: Path,
    ) -> Path:
        """Run Kallisto quantification on individual sequencing library file.

        Args:
            fastq: Path to FASTQ file.
            index: Path to Kallisto index.

        Returns:
            Path to output directory.

        Raises:
            KallistoProblem: Kallisto quantification failed.
        """
        LOGGER.debug(f"Running Kallisto quantification for: {fastq}")

        with tempfile.TemporaryDirectory(
            prefix="kallisto_",
            dir=self.tmp_dir,
        ) as tmp_dir:
            results_dir = Path(tmp_dir)

        cmd = [
            "kallisto",
            "quant",
            "--single",
            "--fragment-length", str(100),
            "--sd", str(300),
            "--index", f"{str(index)}",
            "--output-dir", f"{str(results_dir)}",
            f"{str(fastq)}",
        ]

        result = sp.run(
            cmd,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            LOGGER.error(result.stderr)
            raise KallistoProblem("Failed to run Kallisto index")

        LOGGER.debug(f"Kallisto quantification available at: {results_dir}")
        return results_dir

    @staticmethod
    def get_source_expression(
        kallisto_dir: Path,
    ) -> Dict[Tuple[str, int], float]:
        """Return percentages of total expression per read source.

        Args:
            kallisto_dir: Directory containing Kallisto quantification results.

        Returns:
            Dictionary with percentages of total expression per read source.
                Keys are tuples of source short names and taxon identifiers.

        Raises:
            FileProblem: Kallisto quantification results could not be
                processed.
        """
        tpm: DefaultDict[Tuple[str, int], float] = defaultdict(
            lambda: 0.0
        )
        total_tpm: float = 0.0

        # read Kallisto quantification output table
        _file = kallisto_dir / "abundance.tsv"
        try:
            dat = pd.read_csv(
                _file,
                sep='\t',
            )
            for row in range(dat.shape[0]):
                id_val = dat['target_id'][row]  # pylint: disable=E1136
                tpm_val = float(dat['tpm'][row])  # pylint: disable=E1136
                fields = list(
                    map(str, id_val.split("|"))
                )
                short_name = str(fields[3])
                tax_id = int(fields[4])
                tpm[(short_name, tax_id)] += tpm_val
                total_tpm += tpm_val

            # calculating TPM percentages
            if total_tpm != 0:
                for val in tpm.values():
                    val = round((val / total_tpm) * 100, 2)

        except OSError as exc:
            raise FileProblem(
                "Could not process file: {_file}"
            ) from exc

        tpm_dict = dict(tpm)
        return tpm_dict
