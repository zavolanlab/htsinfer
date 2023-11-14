"""Infer library source from sample data."""

import logging
from pathlib import Path
import subprocess as sp
import tempfile

from Bio import SeqIO  # type: ignore
import pandas as pd  # type: ignore
from pandas import DataFrame  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    KallistoProblem,
    TranscriptsFastaProblem,
    UnsupportedSampleSourceException,
)
from htsinfer.models import (
    ResultsSource,
    Source,
    Config,
)
from htsinfer.utils import (
    validate_top_score,
)


LOGGER = logging.getLogger(__name__)


class GetLibSource:
    """Determine the source of FASTQ sequencing of a single- or paired-end
    seguencing library.

    Args:
        config: Container class for all arguments used in inference
                and results produced by the class.

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
        tax_id: Taxonomy ID of the organism.
    """
    def __init__(  # pylint: disable=E1101
        self,
        config: Config,
    ):
        """Class contructor."""
        self.paths = (config.args.path_1_processed,
                      config.args.path_2_processed)
        self.transcripts_file = config.args.t_file_processed
        self.out_dir = config.args.out_dir
        self.tmp_dir = config.args.tmp_dir
        self.min_match_pct = config.args.lib_source_min_match_pct
        self.min_freq_ratio = config.args.lib_source_min_freq_ratio
        self.tax_id = config.args.tax_id

    def evaluate(self) -> ResultsSource:
        """Infer read source.

        Returns:
            Source results object.
        """
        source = ResultsSource()
        # Check if library_source is provided, otherwise infer it
        if self.tax_id is not None:
            source.file_1.taxon_id = self.tax_id
            src_name = self.get_source_name(
                self.tax_id,
                self.transcripts_file
            )
            source.file_1.short_name = src_name

            if self.paths[1] is not None:
                source.file_2.taxon_id = self.tax_id
                source.file_2.short_name = source.file_1.short_name

        else:
            index = self.create_kallisto_index()
            library_source = self.get_source(
                fastq=self.paths[0],
                index=index,
            )
            source.file_1.short_name = library_source.short_name
            source.file_1.taxon_id = library_source.taxon_id

            if self.paths[1] is not None:
                library_source = self.get_source(
                    fastq=self.paths[1],
                    index=index,
                )
                source.file_2.short_name = library_source.short_name
                source.file_2.taxon_id = library_source.taxon_id

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
        LOGGER.debug(f"Determining source for library: {fastq}")
        source: Source = Source()

        # run quantification
        kallisto_dir = self.run_kallisto_quantification(
            fastq=fastq,
            index=index,
        )

        # process expression levels
        tpm_df = self.get_source_expression(
            kallisto_dir=kallisto_dir,
        )

        # write data frame (in JSON) to file
        filename = (
            Path(self.out_dir) / f"library_source_{fastq.name}.json"
        )
        LOGGER.debug(f"Writing results to file: {filename}")
        tpm_df.to_json(
            filename,
            orient='split',
            index=False,
            indent=True,
        )

        # validate results
        if validate_top_score(
            vector=tpm_df['tpm'].to_list(),
            min_value=self.min_match_pct,
            min_ratio=self.min_freq_ratio,
            rev_sorted=True,
            accept_zero=True,
        ):
            source.short_name, taxon_id = tpm_df.iloc[0]['source_ids']
            source.taxon_id = int(taxon_id)

        LOGGER.debug(f"Source: {source}")
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
            if "zero reads pseudoaligned" not in str(result.stderr):
                LOGGER.error(result.stderr)
                raise KallistoProblem("Failed to run Kallisto quantification")

        LOGGER.debug(f"Kallisto quantification available at: {results_dir}")
        return results_dir

    @staticmethod
    def get_source_expression(
        kallisto_dir: Path,
    ) -> DataFrame:
        """Return percentages of total expression per read source.

        Args:
            kallisto_dir: Directory containing Kallisto quantification results.

        Returns:
            Data frame with columns `source_ids` (a tuple of source short name
                and taxon identifier, e.g., `("hsapiens", 9606)`) and `tpm`,
                signifying the percentages of total expression per read source.
                The data frame is sorted by total expression in descending
                order.

        Raises:
            FileProblem: Kallisto quantification results could not be
                processed.
        """
        # pylint: disable=E1101,E1137

        # read Kallisto quantification output table
        _file = kallisto_dir / "abundance.tsv"
        try:
            dat = pd.read_csv(
                _file,
                sep='\t',
            )
        except OSError as exc:
            raise FileProblem(
                "Could not process file: {_file}"
            ) from exc

        if dat.empty:
            raise TranscriptsFastaProblem(
                "Empty abundance.tsv file created by kallisto quantification."
            )

        # handle case where no alignments are found
        dat.tpm.fillna(0, inplace=True)

        # aggregate expression by source identifiers
        dat[[
            'gene_symbol',
            'gene_id',
            'transcript_id',
            'short_name',
            'taxon_id'
        ]] = dat.target_id.str.split('|', n=4, expand=True)
        dat['source_ids'] = list(zip(dat.short_name, dat.taxon_id))
        total_tpm = dat.tpm.sum()
        dat_agg = dat.groupby(['source_ids'])[['tpm']].agg('sum')
        dat_agg['source_ids'] = dat_agg.index
        dat_agg.reset_index(drop=True, inplace=True)

        # calculate percentages
        if total_tpm != 0:
            dat_agg.tpm = dat_agg.tpm / total_tpm * 100

        # return as dictionary
        return dat_agg.sort_values(["tpm"], ascending=False)

    @staticmethod
    def get_source_name(
        taxon_id: int,
        transcripts_file: Path,
    ) -> str:
        """Return name of the source organism, based on tax ID.

        Args:
            taxon_id: Taxonomy ID of a given organism.
            transcripts_file: Path to FASTA file containing transcripts.

        Returns:
            Short name of the organism belonging to the given tax ID.

        Raises:
            FileProblem: Could not process input FASTA file.
            UnsupportedSampleSourceException: Taxon ID is not supported.
        """
        src_dict = {}

        try:
            for record in list(SeqIO.parse(
                    handle=transcripts_file,
                    format='fasta',
            )):
                tax_id = int(record.description.split("|")[4])
                src_name = record.description.split("|")[3]

                src_dict[tax_id] = src_name

        except OSError as exc:
            raise FileProblem(
                f"Could not process file '{transcripts_file}'"
            ) from exc

        try:
            return src_dict[taxon_id]

        except KeyError as exc:
            raise UnsupportedSampleSourceException(
                f'Taxon ID "{taxon_id}" is not supported by HTSinfer.'
            ) from exc
