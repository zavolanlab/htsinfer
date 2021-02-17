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
    max_records: int = 100000
) -> None:
    """Infers adapter information for one or two fastq files.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record.
     """
    # Process file 1
    logger.debug(f"Processing file 1: {file_1}")
    process_fastq_file(file_1, max_records)

    # Process file 2
    if file_2:
        logger.debug(f"Processing file 2: {file_2}")
        process_fastq_file(file_2, max_records)


def process_fastq_file(
    file: str,
    max_records: int = 10000
) -> None:
    """Process adapters count info.

    Args:
        file (str) : File path to read/first mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record.
    """
    _open = partial(
        gzip.open, mode='rt'
    ) if file.endswith(".gz") else open

    try:
        logger.debug("Opening file...")
        with _open(file) as _file:

            adapter_counts: Dict[str, int] = {}
            records: int = 0
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

                # Update records
                records += 1

                # End if requested records processed
                if max_records and records >= max_records:
                    break

            # Converting dictionary into dataframe
            adapters_df = pd.DataFrame(adapter_counts.items())
            adapters_df.columns = ['Adapter', 'Count']
            adapters_df = adapters_df.sort_values(by='Count', ascending=False).reset_index(drop=True)
            adapters_df.to_csv(f"{file}_adapters_count.csv")

    except ValueError as exc:
        logger.error(f"Invalid input file '{file}'. Error: {str(exc)}")
