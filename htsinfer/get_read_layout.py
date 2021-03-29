"""Infer adapter sequences present in reads."""

from collections import defaultdict
import logging
from pathlib import Path
from typing import (DefaultDict, List, Optional, Tuple)

import ahocorasick as ahc  # type: ignore
from Bio.SeqIO.QualityIO import FastqGeneralIterator  # type: ignore
from pandas import DataFrame  # type: ignore

from htsinfer.exceptions import FileProblem
from htsinfer.models import ResultsLayout
from htsinfer.utils import minmatch_factor_validator

LOGGER = logging.getLogger(__name__)


class GetReadLayout:
    """Determine the adapter sequence present in the FASTQ sequencing libraries.

    Args:
        adapter_file: Adapter file containing the list of all adapter sequences
            that neeeds to be searched in the FASTQ files.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter sequence to be considered as the
            resulting sequence.
        factor: The minimum frequency ratio between the first and second most
            frequent adapter in order for an adapter sequence to be returned
            as the resulting sequence.
        out_dir: Path to directory where output is written to.

    Attrubutes:
        adapter_file: Adapter file containing the list of all adapter sequences
            that neeeds to be searched in the FASTQ files.
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter sequence to be considered as the
            resulting sequence.
        factor: The minimum frequency ratio between the first and second most
            frequent adapter in order for an adapter sequence to be returned
            as the resulting sequence.
        out_dir: Path to directory where output is written to.
        results: Results container for storing adapter sequence information for
            the provided files.
    """
    def __init__(
        self,
        adapter_file: Path,
        path_1: Path,
        path_2: Optional[Path] = None,
        min_match: float = 5,
        factor: float = 2,
        out_dir: Path = Path.cwd(),
    ):
        """Class contructor."""
        self.adapter_file: Path = adapter_file
        self.path_1: Path = path_1
        self.path_2: Optional[Path] = path_2
        self.min_match: float = min_match
        self.factor: float = factor
        self.out_dir: Path = out_dir
        self.results: ResultsLayout = ResultsLayout()

    def evaluate(self) -> None:
        """Decide adapter sequence."""

        # process file 1
        LOGGER.debug(f"Processing file: '{self.path_1}'")
        adapter_file_1 = GetAdapter(
            adapter_file=self.adapter_file,
            path=self.path_1,
            min_match=self.min_match,
            factor=self.factor,
            out_dir=self.out_dir,
        )
        adapter_file_1.evaluate()
        self.results.file_1 = adapter_file_1.result
        LOGGER.debug(f"Adapter Sequence: {self.results.file_1}")

        # process file 2
        if self.path_2 is not None:
            LOGGER.debug(f"Processing file: '{self.path_2}'")
            adapter_file_2 = GetAdapter(
                adapter_file=self.adapter_file,
                path=self.path_2,
                min_match=self.min_match,
                factor=self.factor,
                out_dir=self.out_dir,
            )
            adapter_file_2.evaluate()
            self.results.file_2 = adapter_file_2.result
            LOGGER.debug(f"Adapter Sequence: {self.results.file_2}")


class GetAdapter():
    """Determine adapter sequence for an individual FASTQ library.

    Args:
        adapter_file: Adapter file containing the list of all adapter sequences
            that neeeds to be searched in the FASTQ files.
        path: File path to read library.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter sequence to be considered as the
            resulting sequence.
        factor: The minimum frequency ratio between the first and second most
            frequent adapter in order for an adapter sequence to be returned
            as the resulting sequence.
        out_dir: Path to directory where output is written to.

    Attributes:
        adapter_file: Adapter file containing the list of all adapter sequences
            that neeeds to be searched in the FASTQ files.
        path: File path to read library.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter sequence to be considered as the
            resulting sequence.
        factor: The minimum frequency ratio between the first and second most
            frequent adapter in order for an adapter sequence to be returned
            as the resulting sequence.
        out_dir: Path to directory where output is written to.
        adapters: List of adapters sequence.
        trie: Trie datastructure of adapters sequence.
        adapter_counts: Dictionary of adapter sequence and its count
            percentage.
        result: The most frequent adapter sequence in FASTQ file.
    """
    def __init__(
        self,
        adapter_file: Path,
        path: Path,
        min_match: float = 5,
        factor: float = 2,
        out_dir: Path = Path.cwd(),
    ):
        """Class constructor."""
        self.adapter_file: Path = adapter_file
        self.path: Path = path
        self.min_match: float = min_match
        self.factor: float = factor
        self.out_dir: Path = out_dir
        self.adapters: List[Tuple[str, int]] = []
        self.trie: ahc.Automaton = ahc.Automaton()
        self.adapter_counts: DefaultDict[str, float] = defaultdict(lambda: 0)
        self.result: Optional[str] = None

    def evaluate(self) -> None:
        """Search for adapters sequence and validates minimum match percentage
        and minimum frequency ratio."""
        # load adapters
        self._load_adapters()

        # creates trie of adapters
        self._make_aho_auto()

        records: int = 0
        total_count: int = 0

        try:
            with open(self.path) as _f:  # type: ignore

                LOGGER.debug("Procecssing Reads")
                try:
                    for record in FastqGeneralIterator(source=_f):
                        read = record[1]

                        # Searching for adapters in read
                        for _, (_, key) in self.trie.iter(read):
                            self.adapter_counts[key] += 1
                            total_count += 1

                        records += 1

                except StopIteration as exc:
                    self.result = None
                    raise FileProblem(f"File is corrupt: {self.path}") from exc
        except (OSError, ValueError) as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        LOGGER.debug(f"Total records processed: {records}")

        if total_count != 0:
            for i in self.adapter_counts:
                self.adapter_counts[i] = round(
                    (self.adapter_counts[i]/records)*100, 2
                    )
            adapters_df = self._covert_dic_to_df()

            # Checks validator
            if minmatch_factor_validator(
                adapters_df,
                column_index=1,
                min_match=self.min_match,
                factor=self.factor,
            ):
                self.result = adapters_df.iloc[0]['Adapter Sequence']
            else:
                self.result = None
        else:
            self.result = None

    def _covert_dic_to_df(self) -> DataFrame:
        """Converting dictionary into dataframe and writing json file.

        Returns:
            Dataframe of adapter sequence with and its count percentage.
        """
        adapters_df = DataFrame(self.adapter_counts.items())
        adapters_df.columns = ['Adapter Sequence', 'Count Percentage']
        adapters_df = adapters_df.sort_values(
            by='Count Percentage', ascending=False
            ).reset_index(drop=True)
        name = Path(self.out_dir) / f"ReadLayout_{Path(self.path).name}.json"
        LOGGER.debug(f"Creating {name}")
        adapters_df.to_json(
            name, orient='split', index=False, indent=True
            )
        return adapters_df

    def _make_aho_auto(self) -> None:
        """Adding all adapters sequence into trie datastructure."""
        LOGGER.debug("Creating trie of adapters")
        for (adapter, tag) in self.adapters:
            self.trie.add_word(adapter, (tag, adapter))

        self.trie.make_automaton()

    def _load_adapters(self) -> None:
        """Loads adapters sequence from the adapters_list.txt file."""
        LOGGER.debug("Loading adapters")
        _file = open(self.adapter_file, "r")
        tag: int = 1
        while True:
            line = _file.readline()
            name = line.strip()
            if not name:
                break
            self.adapters.append((name, tag))
            tag += 1
