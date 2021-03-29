"""Infer organism information for sequencing library."""

import logging
from pathlib import Path
import subprocess as sp
import tempfile
from typing import (Dict, Optional, Tuple)
import zipfile as zp

import pandas as pd  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    KallistoProblem,
)
from htsinfer.models import ResultsSource
from htsinfer.utils import minmatch_factor_validator

LOGGER = logging.getLogger(__name__)


class GetLibSource:
    """Determine the organism present in the FASTQ sequencing libraries.
    Args:
        fasta: File path to transcripts FASTA file.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        min_match: Minimum percentage that given organism needs to have
            to be considered as the resulting organism.
        factor: The minimum frequency ratio between the first and second
            most frequent organism in order for organism to be
            considered as the resulting organism.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.

    Attrubutes:
        fasta: File path to transcripts FASTA file.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        min_match: Minimum percentage that given organism needs to have
            to be considered as the resulting organism.
        factor: The minimum frequency ratio between the first and second
            most frequent organism in order for organism to be
            considered as the resulting organism.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        results: Results container for storing organism information for
            the provided files.
    """
    def __init__(
        self,
        fasta: Path,
        path_1: Path,
        path_2: Optional[Path] = None,
        min_match: float = 5,
        factor: float = 2,
        out_dir: Path = Path.cwd(),
        tmp_dir: Path = Path(tempfile.gettempdir()),
    ):
        """Class contructor."""
        self.fasta: Path = fasta
        self.path_1: Path = path_1
        self.path_2: Optional[Path] = path_2
        self.min_match: float = min_match
        self.factor: float = factor
        self.out_dir: Path = out_dir
        self.tmp_dir = tmp_dir
        self.results: ResultsSource = ResultsSource()

    def evaluate(self) -> None:
        """Decides organism."""
        # Extracts transcripts fasta file
        try:
            with zp.ZipFile(self.fasta, "r") as zip_ref:
                zip_ref.extractall(self.tmp_dir)
        except (FileNotFoundError, zp.BadZipFile) as exc:
            self.results.file_1 = None
            self.results.file_2 = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        # Runs Kallisto index
        try:
            self._kallisto_index()
        except KallistoProblem as exc:
            self.results.file_1 = None
            self.results.file_2 = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        # process file 1
        LOGGER.debug(f"Processing file: '{self.path_1}'")
        organism_file_1 = GetOrganism(
            path=self.path_1,
            min_match=self.min_match,
            factor=self.factor,
            out_dir=self.out_dir,
            tmp_dir=self.tmp_dir,
        )
        organism_file_1.evaluate()
        self.results.file_1 = organism_file_1.result
        LOGGER.debug(f"Organism: {self.results.file_1}")

        # process file 2
        if self.path_2 is not None:
            LOGGER.debug(f"Processing file: '{self.path_2}'")
            organism_file_2 = GetOrganism(
                path=self.path_2,
                min_match=self.min_match,
                factor=self.factor,
                out_dir=self.out_dir,
                tmp_dir=self.tmp_dir,
            )
            organism_file_2.evaluate()
            self.results.file_2 = organism_file_2.result
            LOGGER.debug(f"Organism: {self.results.file_2}")

    def _kallisto_index(self) -> None:
        """Builds an index from a FASTA formatted file of target sequences."""
        LOGGER.debug("Running Kallisto index")
        _file = Path(Path(self.fasta).name).stem
        index_cmd = "kallisto index -i " + str(self.tmp_dir) + \
            "/transcripts.idx --make-unique " + str(self.tmp_dir) + "/" + _file
        result = sp.run(index_cmd, shell=True, capture_output=True, text=True)
        LOGGER.debug(result.stderr)
        if result.returncode == 0:
            pass
        else:
            raise KallistoProblem("Failed to run Kallisto index")


class GetOrganism():
    """Determine organism for an individual FASTQ library.
    Args:
        path: File path to read library.
        min_match: Minimum percentage that given organism needs to have
            to be considered as the resulting organism.
        factor: The minimum frequency ratio between the first and second
            most frequent organism in order for organism to be
            considered as the resulting organism.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.

    Attributes:
        path: File path to read library.
        min_match: Minimum percentage that given organism needs to have
            to be considered as the resulting organism.
        factor: The minimum frequency ratio between the first and second
            most frequent organism in order for organism to be
            considered as the resulting organism.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        organism_count: Dictionary with count percentage for all organisms.
        result: The most frequent organism in FASTQ file.
    """
    def __init__(
        self,
        path: Path,
        min_match: float = 5,
        factor: float = 2,
        out_dir: Path = Path.cwd(),
        tmp_dir: Path = Path(tempfile.gettempdir()),
    ):
        """Class constructor."""
        self.path: Path = path
        self.min_match: float = min_match
        self.factor: float = factor
        self.out_dir: Path = out_dir
        self.tmp_dir: Path = tmp_dir
        self.organism_count: Dict[Tuple[str, int], float] = {}
        self.result: Optional[str] = None

    def evaluate(self) -> None:
        """Determines organism and validates minimum match percentage
        and minimum frequency ratio."""
        # Runs Kallisto quant
        try:
            self._kallisto_quant()
        except KallistoProblem as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        self._process_count_info()
        # Converts dictionary to dataframe
        organism_df = self._convert_dic_to_df()

        # Checks validator
        if minmatch_factor_validator(
            count_df=organism_df,
            column_index=0,
            min_match=self.min_match,
            factor=self.factor,
        ):
            self.result = organism_df.iloc[0]['Organism']
        else:
            self.result = None

    def _kallisto_quant(self) -> None:
        """Runs the quantification algorithm."""
        LOGGER.debug("Running Kallisto quant")
        quant_cmd = "kallisto quant -i " + str(self.tmp_dir) + \
            "/transcripts.idx -o " + str(self.tmp_dir) + \
            "/output -l 100 -s 300 --single " + str(self.path)
        result = sp.run(quant_cmd, shell=True, capture_output=True, text=True)
        LOGGER.debug(result.stderr)
        if result.returncode == 0:
            pass
        else:
            raise KallistoProblem("Failed to run Kallisto quant")

    def _process_count_info(self) -> None:
        """Process organisms count information."""
        # Reading tsv file created by kallisto quant
        try:
            _file = Path(self.tmp_dir) / "output" / "abundance.tsv"
            abundance_df = pd.read_csv(_file, sep='\t')
            dimension = abundance_df.shape
            rows = dimension[0]
            total_tpm: float = 0.0
            for i in range(rows):
                row_id = abundance_df['target_id'][i]
                contents = list(map(str, row_id.split("|")))
                organism_name = contents[3]
                organism_tax_id = int(contents[4])
                # Update organism tpm count
                if (organism_name, organism_tax_id) in self.organism_count:
                    self.organism_count[(organism_name, organism_tax_id)] += \
                        float(abundance_df['tpm'][i])
                else:
                    self.organism_count[(organism_name, organism_tax_id)] = \
                        float(abundance_df['tpm'][i])
                total_tpm += float(abundance_df['tpm'][i])

            # Calculating Percentage
            if total_tpm != 0:
                for (_name, _id), _ in self.organism_count.items():
                    self.organism_count[(_name, _id)] = round(
                        (self.organism_count[(_name, _id)]/total_tpm)*100, 2
                        )
        except (FileNotFoundError, OSError) as exc:
            self.result = None
            raise KallistoProblem(
                "Could not open abundance.tsv created by Kallisto quant"
                ) from exc

    def _convert_dic_to_df(self) -> pd.DataFrame:
        """Converting dictionary into dataframe and writing json file.

        Returns:
            Dataframe of count info percentage for all organisms.
        """
        organism_df = pd.DataFrame(self.organism_count.items())
        organism_df[['Organism', 'Taxon ID']] = pd.DataFrame(
            organism_df[0].tolist()
            )
        organism_df = organism_df.sort_values(
            by=1, ascending=False
            ).reset_index(drop=True).drop([0], axis=1)
        organism_df = organism_df.rename(columns={1: 'Match Percentage'})
        name = Path(self.out_dir) / f"ReadSource_{Path(self.path).name}.json"
        LOGGER.debug(f"Creating {name}")
        organism_df.to_json(
            name, orient='split', index=False, indent=True
            )
        return organism_df
