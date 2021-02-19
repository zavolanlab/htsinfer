"""Infer adapter sequences present in reads."""

from functools import partial
import pandas as pd
import ahocorasick as ahc
import logging
import gzip
from typing import (Dict, List, Tuple)

from Bio.SeqIO.QualityIO import FastqGeneralIterator

logger = logging.getLogger(__name__)


adapters = [
    ('AAAAAAAAAAAAAAA', 1),
    ('AATGATACGGCGACC', 2),
    ('ACACTCTTTCCCTAC', 3),
    ('ACAGGTTCAGAGTTC', 4),
    ('ATCTCGTATGCCGTC', 5),
    ('CAAGCAGAAGACGGC', 6),
    ('CCGACAGGTTCAGAG', 7),
    ('CCGAGCCCACGAGAC', 8),
    ('CCGGTTCTTCCCTGC', 9),
    ('CGACAGGTTCAGAGT', 10),
    ('CGGTCTCGGCATTCC', 11),
    ('CTGTAGGCACCATCA', 12),
    ('GAAUUCCACCACGUU', 13),
    ('GATCGGAAGAGCACA', 14),
    ('GATCGGAAGAGCGGT', 15),
    ('GATCGGAAGAGCTCG', 16),
    ('GATCGTCGGACTGTA', 17),
    ('GCCTTGGCACCCGAG', 18),
    ('GTCTCGTGGGCTCGG', 19),
    ('GTGACTGGAGTTCAG', 20),
    ('GUUCAGAGUUCUACA', 21),
    ('TCGGACTGTAGAACT', 22),
    ('TCGTATGCCGTCTTC', 23),
    ('TCGTCGGCAGCGTCA', 24),
    ('TGGAATTCTCGGGTG', 25),
    ('UCGUAUGCCGUCUUC', 26),
    ('AGATCGGAAGAGCGT', 27)
]


def make_aho_auto(
    adapters: List[Tuple[str, int]]
):
    logger.debug("Creating trie")
    trie = ahc.Automaton()
    for (adapter, tag) in adapters:
        trie.add_word(adapter, (tag, adapter))

    trie.make_automaton()
    return trie


trie = make_aho_auto(adapters)


def infer(
    file_1: str,
    file_2: str = None,
    max_records: int = 100000,
    min_match: float = 10,
    factor: float = 2
) -> Tuple[str, str]:
    """Infers adapter information for one or two fastq files.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record.
        min_match (float) : minimum match percentage that first adapter needs
        to have.
        factor (float) : factor by which first adapter is greater than the
        second adapter.

    Returns:
        Type of Adapter that is present in library.
    """
    # Process file 1
    logger.debug(f"Processing file 1: {file_1}")
    result_1 = process_fastq_file(file_1, max_records, min_match, factor)

    # Process file 2
    result_2 = "not_available"
    if file_2:
        logger.debug(f"Processing file 2: {file_2}")
        result_2 = process_fastq_file(file_2, max_records, min_match, factor)

    logger.debug("Returning results...")
    return (result_1, result_2)


def process_fastq_file(
    file: str,
    max_records: int = 10000,
    min_match: float = 10,
    factor: float = 2
) -> str:
    """Process adapters count info.

    Args:
        file (str) : File path to read/first mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record.
        min_match (float) : minimum match percentage that first adapter needs
        to have.
        factor (float) : factor by which first adapter is greater than the
        second adapter.

    Returns:
        Adapter information
    """
    _open = partial(
        gzip.open, mode='rt'
    ) if file.endswith(".gz") else open

    try:
        logger.debug("Opening file...")
        with _open(file) as _file:

            adapter_counts: Dict[str, int] = {}
            records: int = 0
            total_count = 0
            for record in FastqGeneralIterator(source=_file):
                # Get next read
                read = record[1]
                logger.debug("Processing read")
                # Searching for adapters in read
                for loc, (tag, key) in trie.iter(read):
                    if key in adapter_counts:
                        adapter_counts[key] += 1
                    else:
                        adapter_counts[key] = 1
                    total_count += 1

                # Update records
                records += 1

                # End if requested records processed
                if max_records and records >= max_records:
                    break

            if total_count == 0:
                return "NA"

            # Calculating Percentage
            for i in adapter_counts:
                adapter_counts[i] = round((adapter_counts[i]/total_count)*100, 2)

            # Converting dictionary into dataframe
            adapters_df = pd.DataFrame(adapter_counts.items())
            adapters_df.columns = ['Adapter', 'Count %']
            adapters_df = adapters_df.sort_values(by='Count %', ascending=False).reset_index(drop=True)
            logger.debug(f"Creating {file}_adapters_count.csv")
            adapters_df.to_csv(f"{file}_adapters_count.csv")

            # Checking confidence score
            if(confidence(adapters_df, min_match, factor)):
                result = adapters_df.iloc[0]['Adapter']
            else:
                result = "NA"
            return result

    except OSError:
        logger.error(f"Invalid input file '{file}'")
        return "invalid_file"


def confidence(
    adapters_df: pd.DataFrame,
    min_match: float = 10,
    factor: float = 2
) -> bool:
    """Checks confidence score

    Args:
        adapters_df (dataframe) : adapters count.
        min_match (float) : minimum match percentage that first adapter needs
        to have.
        factor (float) : factor by which first adapter is greater than the
        second adapter.

    Returns:
        bool : whether it satisfies confidence score or not.
    """
    if adapters_df.iloc[0]['Count %'] < min_match:
        return False
    if adapters_df.iloc[1]['Count %'] != 0:
        ratio = adapters_df.iloc[0]['Count %']/adapters_df.iloc[1]['Count %']
        if ratio < factor:
            return False
    return True
