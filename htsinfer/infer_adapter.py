"""Infer adapter sequences present in reads."""

from functools import partial
import gzip
import logging
from typing import (Dict, List, Tuple)

import ahocorasick as ahc  # type: ignore
from Bio.SeqIO.QualityIO import FastqGeneralIterator  # type: ignore
from pandas import DataFrame  # type: ignore

LOGGER = logging.getLogger(__name__)


def infer(
    adapter_file: str,
    file_1: str,
    file_2: str = None,
    min_match: float = 10,
    factor: float = 2,
    max_records: int = 100000,
) -> Tuple[str, str]:
    """Infers adapter information for one or two fastq files.

    Args:
        adapter_file: Adapter file containing the list of all adapter sequences
            that neeeds to be searched in the FASTQ files.
        file_1: File path to first mate library.
        file_2: File path to second mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter to be considered as the resulting
            adapter. If no adapter is found more frequently than the specified
            number (in percent), null is returned in the JSON result to
            indicate that no adapter could be confidently identified.
        factor: The minimum frequency ratio between the first and second most
            frequent adapters in order for an adapter sequence to be returned
            in the JSON result. If frequency ratio is less than the specified
            value, null is returned in the JSON result to indicate that no
            single adapter could be confidently identified.

    Returns:
        Adapter sequence that is present in library.
    """
    adapters = load_adapters(adapter_file)

    # Reading file 1
    LOGGER.debug(f"Reading mate 1 file: {file_1}")
    result_1 = read_fastq_file(
        adapters, file_1, max_records, min_match, factor
        )

    # Reading file 2
    result_2 = "not_available"
    if file_2 is not None:
        LOGGER.debug(f"Reading mate 2 file: {file_2}")
        result_2 = read_fastq_file(
            adapters, file_2, max_records, min_match, factor
            )

    LOGGER.debug("Returning results...")
    return (result_1, result_2)


def read_fastq_file(
    adapters: List[Tuple[str, int]],
    file: str,
    max_records: int = 10000,
    min_match: float = 10,
    factor: float = 2,
) -> str:
    """Process adapters count info.

    Args:
        adapters: List of adapters sequences.
        file: File path to mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter to be considered as the resulting
            adapter.
        factor: The minimum frequency ratio between the first and second most
            frequent adapters in order for an adapter sequence to be returned
            in the JSON result.

    Returns:
        Adapter sequence
    """
    _open = partial(
        gzip.open, mode='rt'
    ) if file.endswith(".gz") else open

    trie = make_aho_auto(adapters)
    adapter_counts: Dict[str, float] = {}
    records: int = 0
    total_count: int = 0

    try:
        LOGGER.debug("Opening file...")
        with _open(file) as _file:  # type: ignore

            LOGGER.debug("Processing reads")
            for record in FastqGeneralIterator(source=_file):
                read = record[1]

                # Searching for adapters in read
                for _, (_, key) in trie.iter(read):
                    if key in adapter_counts:
                        adapter_counts[key] += 1
                    else:
                        adapter_counts[key] = 1
                    total_count += 1

                records += 1

                # End if requested records processed
                if max_records and records >= max_records:
                    break

    except OSError:
        LOGGER.error(f"Invalid input file '{file}'")
        return "invalid_file"

    if total_count == 0:
        return "NA"

    # Calculating Percentage
    for i in adapter_counts:
        adapter_counts[i] = round(
            (adapter_counts[i]/records)*100, 2
            )

    adapters_df = covert_dic_to_df(adapter_counts, file)
    # Checking confidence score
    if confidence(adapters_df, min_match, factor):
        result = adapters_df.iloc[0]['Adapter']
    else:
        result = "NA"
    return result


def covert_dic_to_df(
    adapter_counts: Dict[str, float],
    file: str
) -> DataFrame:
    """Converting dictionary into dataframe and writing json file.

    Args:
        adapter_counts: Dictionary of adapter sequence and its count
        percentage.
        file: File path to mate library.

    Returns:
        Dataframe of adapter sequence with and its count percentage.

    Examples:
        Adapter                Count %
        GATCGGAAGAGCACA        2.83
        AAAAAAAAAAAAAAA        1.79
    """
    adapters_df = DataFrame(adapter_counts.items())
    adapters_df.columns = ['Adapter', 'Count %']
    adapters_df = adapters_df.sort_values(
        by='Count %', ascending=False
        ).reset_index(drop=True)
    LOGGER.debug(f"Creating {file}_adapters_count.json")
    adapters_df.to_json(f'{file}_adapters_count.json', orient='records')
    return adapters_df


def load_adapters(
    adapter_file: str
) -> List[Tuple[str, int]]:
    """Loads adapters sequence from the adapter.txt file.

    Args:
        adapter_file: Adapter file containing the list of all adapter sequence
            that neeeds to be searched in the FASTQ files.

    Returns:
        List of adapters sequence.
    """

    _file = open(adapter_file, "r")
    adapters: List[Tuple[str, int]] = []
    tag: int = 1
    while True:
        line = _file.readline()
        name = line.strip()
        if not name:
            break
        adapters.append((name, tag))
        tag += 1

    return adapters


def make_aho_auto(
    adapters: List[Tuple[str, int]]
) -> ahc.Automaton:
    """Adding all adapters sequence into trie datastructure.

    Args:
        adapters: List of adapters sequence.

    Returns:
        Returns trie of adapters sequence.
    """
    LOGGER.debug("Creating trie")
    trie = ahc.Automaton()
    for (adapter, tag) in adapters:
        trie.add_word(adapter, (tag, adapter))

    trie.make_automaton()
    return trie


def confidence(
    adapters_df: DataFrame,
    min_match: float = 10,
    factor: float = 2
) -> bool:
    """Checks confidence score

    Args:
        adapters_df: Adapters sequence count.
        min_match: Minimum percentage of reads that contain a given adapter
            in order for that adapter to be considered as the resulting
            adapter.
        factor: The minimum frequency ratio between the first and second most
            frequent adapters in order for an adapter sequence to be returned
            in the JSON result.

    Returns:
        Whether it satisfies confidence score or not.
    """
    if adapters_df.iloc[0]['Count %'] < min_match:
        return False
    if adapters_df.iloc[1]['Count %'] != 0:
        ratio = adapters_df.iloc[0]['Count %']/adapters_df.iloc[1]['Count %']
        if ratio < factor:
            return False
    return True
