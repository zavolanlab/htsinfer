"""Infer adapter sequences present in reads."""

import logging
from pathlib import Path
from typing import (Dict, List, Optional, Tuple)

import ahocorasick as ahc  # type: ignore
from Bio.SeqIO.QualityIO import FastqGeneralIterator  # type: ignore

from htsinfer.exceptions import FileProblem
from htsinfer.models import ResultsLayout
from htsinfer.utils import (
    convert_dict_to_df,
    validate_top_score,
)

LOGGER = logging.getLogger(__name__)


class GetReadLayout:
    """Determine the adapter sequence present in the FASTQ sequencing libraries.

    Args:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        adapter_file: Path to text file containing 3' adapter sequences (one
            sequence per line) to scan for.
        out_dir: Path to directory where output is written to.
        min_match_pct: Minimum percentage of reads that contain a given adapter
            sequence in order for it to be considered as the library's 3'-end
            adapter.
        min_freq_ratio: Minimum frequency ratio between the first and second
            most frequent adapter in order for the former to be considered as
            the library's 3'-end adapter.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        adapter_file: Path to text file containing 3' adapter sequences (one
            sequence per line) to scan for.
        out_dir: Path to directory where output is written to.
        min_match_pct: Minimum percentage of reads that contain a given adapter
            sequence in order for it to be considered as the library's 3'-end
            adapter.
        min_freq_ratio: Minimum frequency ratio between the first and second
            most frequent adapter in order for the former to be considered as
            the library's 3'-end adapter.
        results: Results container for storing adapter sequence information for
            the provided files.

    Examples:
        >>> GetReadLayout(
        ...     path_1="tests/files/sra_sample_2.fastq",
        ...     adapter_file="data/adapter_fragments.txt",
        ... ).evaluate()
        ResultsLayout(
            file_1=<Layout().adapt_3: "AAAAAAAAAAAAAAA">,
            file_2=<Layout().adapt_3: None>,
        )
        >>> GetReadLayout(
        ...     path_1="tests/files/sra_sample_1.fastq",
        ...     path_2="tests/files/sra_sample_2.fastq",
        ...     adapter_file="data/adapter_fragments.txt",
        ...     min_match_pct=2,
        ...     min_freq_ratio=1,
        ... ).evaluate()
        ResultsLayout(
            file_1=<Layout().adapt_3: "AAAAAAAAAAAAAAA">,
            file_2=<Layout().adapt_3: "AAAAAAAAAAAAAAA">,
        )
    """
    def __init__(
        self,
        path_1: Path,
        path_2: Optional[Path] = None,
        adapter_file: Path = (
            Path(__file__).parent.parent.absolute() /
            "data/adapter_fragments.txt"
        ),
        out_dir: Path = Path.cwd() / 'results_htsinfer',
        min_match_pct: float = 2,
        min_freq_ratio: float = 2,
    ):
        """Class contructor."""
        self.path_1: Path = path_1
        self.path_2: Optional[Path] = path_2
        self.adapter_file: Path = adapter_file
        self.out_dir: Path = out_dir
        self.min_match_pct: float = min_match_pct
        self.min_freq_ratio: float = min_freq_ratio
        self.results: ResultsLayout = ResultsLayout()

    def evaluate(self) -> None:
        """Decide adapter sequence."""

        # process file 1
        LOGGER.debug(f"Processing file: '{self.path_1}'")
        file_1_adapt_3 = GetAdapter3(
            path=self.path_1,
            adapter_file=self.adapter_file,
            out_dir=self.out_dir,
            min_match_pct=self.min_match_pct,
            min_freq_ratio=self.min_freq_ratio,
        )
        file_1_adapt_3.evaluate()
        self.results.file_1.adapt_3 = file_1_adapt_3.result
        LOGGER.debug(f"3' adapter sequence: {self.results.file_1}")

        # process file 2
        if self.path_2 is not None:
            LOGGER.debug(f"Processing file: '{self.path_2}'")
            file_2_adapt_3 = GetAdapter3(
                path=self.path_2,
                adapter_file=self.adapter_file,
                out_dir=self.out_dir,
                min_match_pct=self.min_match_pct,
                min_freq_ratio=self.min_freq_ratio,
            )
            file_2_adapt_3.evaluate()
            self.results.file_2.adapt_3 = file_2_adapt_3.result
            LOGGER.debug(f"3' adapter sequence: {self.results.file_2}")


class GetAdapter3():
    """Determine 3' adapter sequence for an individual FASTQ library.

    Args:
        path: File path to read library.
        adapter_file: Path to text file containing 3' adapter sequences (one
            sequence per line) to scan for.
        out_dir: Path to directory where output is written to.
        min_match_pct: Minimum percentage of reads that contain a given adapter
            Minimum percentage of reads that contain a given adapter sequence
            in order for it to be considered as the library's 3'-end adapter.
        min_freq_ratio: Minimum frequency ratio between the first and second
            most frequent adapter in order for the former to be considered as
            the library's 3'-end adapter.

    Attributes:
        path: File path to read library.
        adapter_file: Path to text file containing 3' adapter sequences (one
            sequence per line) to scan for.
        out_dir: Path to directory where output is written to.
        min_match_pct: Minimum percentage of reads that contain a given adapter
            Minimum percentage of reads that contain a given adapter sequence
            in order for it to be considered as the library's 3'-end adapter.
        min_freq_ratio: Minimum frequency ratio between the first and second
            most frequent adapter in order for the former to be considered as
            the library's 3'-end adapter.
        adapters: List of adapter sequences.
        trie: Trie data structure of adapter sequences.
        adapter_counts: Dictionary of adapter sequences and corresponding count
            percentages.
        result: The most frequent adapter sequence in FASTQ file.

    Examples:
        >>> GetAdapter3(
        ...     path_1="tests/files/sra_sample_2.fastq",
        ...     adapter_file="data/adapter_fragments.txt",
        ... ).evaluate()
        <"AAAAAAAAAAAAAAA">
    """
    def __init__(
        self,
        path: Path,
        adapter_file: Path,
        out_dir: Path = Path.cwd(),
        min_match_pct: float = 2,
        min_freq_ratio: float = 2,
    ):
        """Class constructor."""
        self.path: Path = path
        self.adapter_file: Path = adapter_file
        self.out_dir: Path = out_dir
        self.min_match_pct: float = min_match_pct
        self.min_freq_ratio: float = min_freq_ratio
        self.adapters: List[Tuple[str, int]] = []
        self.trie: ahc.Automaton = ahc.Automaton()
        self.adapter_counts: Dict[str, float] = {}
        self.result: Optional[str] = None

    def evaluate(self) -> None:
        """Search for adapter sequences and validate result confidence
        constraints.
        """
        # load adapters
        try:
            self._load_adapters()
        except OSError as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        # create trie of adapters
        self._make_aho_auto()

        # set default counts per adapter
        self.adapter_counts = dict.fromkeys(
            [tup[0] for tup in self.adapters], 0
        )

        records: int = 0
        total_count: int = 0

        try:
            with open(self.path, encoding="utf-8") as _f:  # type: ignore

                LOGGER.debug("Procecssing Reads")
                try:
                    for record in FastqGeneralIterator(source=_f):
                        read = record[1]

                        # Searching for adapters in read
                        for _, (_, key) in self.trie.iter(read):
                            self.adapter_counts[key] += 1
                            total_count += 1

                        records += 1

                except Exception as exc:
                    self.result = None
                    raise FileProblem(f"File is corrupt: {self.path}") from exc

        except (OSError, ValueError) as exc:
            self.result = None
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        LOGGER.debug(f"Total records processed: {records}")

        # process and convert data to data frame
        LOGGER.debug(f"Adapter counts: {dict(self.adapter_counts)}")
        for i in self.adapter_counts:
            self.adapter_counts[i] = round(
                (self.adapter_counts[i]/records)*100, 2
                )
        LOGGER.debug(f"Adapter count fractions: {dict(self.adapter_counts)}")
        adapters_df = convert_dict_to_df(
            dic=self.adapter_counts,
            col_headers=('Adapter Sequence', 'Count Percentage'),
            sort=True,
            sort_by=1,
            sort_ascending=False,
        )

        # write data frame (in JSON) to file
        filename = (
            Path(self.out_dir) / f"read_layout_{Path(self.path).name}.json"
        )
        LOGGER.debug(f"Writing results to file: {filename}")
        adapters_df.to_json(
            filename,
            orient='split',
            index=False,
            indent=True,
        )

        # validate results
        if validate_top_score(
            vector=adapters_df['Count Percentage'].to_list(),
            min_value=self.min_match_pct,
            min_ratio=self.min_freq_ratio,
            rev_sorted=True,
            accept_zero=True,
        ):
            self.result = str(adapters_df.iloc[0]['Adapter Sequence'])
        else:
            self.result = None

    def _make_aho_auto(self) -> None:
        """Add adapter sequences into trie data structure."""
        LOGGER.debug("Creating trie of adapters")
        for (adapter, tag) in self.adapters:
            self.trie.add_word(adapter, (tag, adapter))
        self.trie.make_automaton()

    def _load_adapters(self, allowed_chars="acgtn", max_len=100) -> None:
        """Load adapter sequences from file.

        Args:
            allowed_chars: String composed of accepted/allowed characters in
                each adapter sequence (typically 'acgtn'). It is not necessary
                to provide both lower- and uppercase versions separately.
            max_len: Maximum length of each adapter sequence.

        Raises:
            ValueError: Raised if an adapter sequence is too long or contains
                invalid characters.
        """
        LOGGER.debug("Loading adapters")
        with open(self.adapter_file, encoding="utf-8") as _f:
            i: int = 0
            for line in _f:
                i += 1
                seq = line.strip()
                if not all(
                    c in allowed_chars.lower() for c in seq.lower()
                ):
                    raise ValueError(
                        f"Sequence '{seq}' in adapter file "
                        f"'{self.adapter_file}' contains invalid characters. "
                        f"Only '[{allowed_chars}]' (lower- and uppercase) "
                        "allowed."
                    )
                if len(seq) > max_len:
                    raise ValueError(
                        f"Sequence '{seq}' in adapter file "
                        f"'{self.adapter_file}' is too long.contains. Maximum "
                        f"length is {max_len}."
                    )
                if not seq:
                    continue
                self.adapters.append((seq, i))
