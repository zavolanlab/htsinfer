"""Infer mate information from sample data."""


from enum import Enum
from functools import partial
import gzip
import logging
import re
from typing import (Dict, List, Tuple)

from Bio.SeqIO.QualityIO import FastqGeneralIterator  # type: ignore

LOGGER = logging.getLogger(__name__)


class Outcomes(Enum):
    """Enumerator for mate information outcomes."""
    not_available = "not_available"
    invalid_file = "invalid_file"
    no_mate_info = "no_mate_info"
    first_mate = "first_mate"
    second_mate = "second_mate"
    mixed_mates = "mixed_mates"
    split_mates = "split_mates"
    not_mates = "not_mates"


def infer(
    file_1: str,
    file_2: str = None,
    max_records: int = 10000,
) -> Tuple[str, str, str]:
    """Infers mate information for one or two FASTQ files.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.
        max_records: Limit processing to the indicated number of records,
            starting from the first record. Set to `0` to process all
            records.

    Returns:
        Mate information for each input file and mate relationship.

    Examples:
        >>> infer(file_1="../tests/test_files/first_mate.fastq")
        ('first_mate', 'not_available', 'not_available')
        >>> infer(
        ...     file_1="../tests/test_files/first_mate.fastq",
        ...     file_2="../tests/test_files/second_mate.fastq",
        ...     max_records=10
        ('first_mate', 'second_mate', 'split_mates')
    """
    result_lib_1: str = Outcomes.not_available.value
    result_lib_2: str = Outcomes.not_available.value
    mate_relationship: str = Outcomes.not_available.value

    # Process file 1
    LOGGER.debug(f"Processing file 1: {file_1}")
    result_lib_1, seq_ids_1 = process_fastq_file(
        file=file_1,
        max_records=max_records
    )

    # Process file 2
    if file_2:
        LOGGER.debug(f"Processing file 2: {file_2}")
        result_lib_2, seq_ids_2 = process_fastq_file(
            file=file_2,
            max_records=max_records
        )

        # Check whether libraries are from a pair
        LOGGER.debug("Checking mate relationship between files...")
        if result_lib_1 == "first_mate" and \
                result_lib_2 == "second_mate" and \
                seq_ids_1 == seq_ids_2:
            mate_relationship = Outcomes.split_mates.value
        else:
            mate_relationship = Outcomes.not_mates.value
        LOGGER.debug(f"Mate relationship: {mate_relationship}")

    LOGGER.debug("Returning results...")
    return (result_lib_1, result_lib_2, mate_relationship)


def process_fastq_file(
    file: str,
    max_records: int = 10000,
) -> Tuple[str, set]:
    """Returns mate information for a single FASTQ file.

    Args:
        file: FASTQ file to process.
        max_records: Limit processing to the indicated number of records,
            starting from the first record. Set to `0` to process all
            records.

    Returns:
        Mate information for set of reads and corresponding identifiers.
    """
    seq_ids: Dict[str, List[bool]] = {}

    _open = partial(
        gzip.open, mode='rt'
    ) if file.endswith(".gz") else open

    try:
        LOGGER.debug("Opening file...")
        with _open(file) as _file:  # type: ignore

            records: int = 0
            mate_counts: Dict[int, int] = {1: 0, 2: 0}

            for record in FastqGeneralIterator(source=_file):
                # Get next read
                seq_id = record[0]
                LOGGER.debug(f"Processing read: {seq_id}")

                # Get mate information
                LOGGER.debug("Extracting mate information...")
                mate, seq_id_no_mate = parse_mate_info_from_id(seq_id)
                LOGGER.debug(f"Mate information: {mate}")

                # Check FASTQ consistency
                if seq_id_no_mate not in seq_ids:
                    seq_ids[seq_id_no_mate] = [False, False]
                elif seq_ids[seq_id_no_mate][mate-1]:
                    raise ValueError("Duplicate read names")
                seq_ids[seq_id_no_mate][mate-1] = True

                # Mate information could not be parsed
                if not mate:
                    LOGGER.debug(
                        "Could not extract mate information for read"
                        f"identifier '{seq_id}'."
                    )
                    return (Outcomes.no_mate_info.value, set(seq_ids.keys()))

                # Update counts
                mate_counts[mate] += 1
                records += 1

                # End if requested records processed
                if max_records and records >= max_records:
                    break

            LOGGER.debug(
                "Summarizing mate info..."
            )
            result = summarize_mate_info(counts=mate_counts)
            LOGGER.debug(
                f"Mate info: {result}"
            )

            return (result, set(seq_ids.keys()))

    except OSError:
        LOGGER.error(
            f"Invalid input file '{file}'. Error: Could not open file"
        )
        return (Outcomes.invalid_file.value, set(seq_ids.keys()))

    except ValueError as exc:
        LOGGER.error(f"Invalid input file '{file}'. Error: {str(exc)}")
        return (Outcomes.invalid_file.value, set(seq_ids.keys()))


def parse_mate_info_from_id(
    seq_id: str
) -> Tuple[int, str]:
    """Extracts mate information from sequence identifier.

    Args:
        seq_id: Sequence identifier.

    Returns:
        An integer representing the mate information and the sequence
            identifier minus the mate information, if it could be extracted,
            or the original sequence identifier otherwise.
            The following mate information integers may be observed:
                0: Mate information could not be extracted
                1: First mate
                2: Second mate
    """
    # Store sequence ID conventions and index of correspoding capture group
    conventions: List[re.Pattern] = []
    # Illumina Casava >=1.8
    # example EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    # https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    # https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
    conventions.append(
        re.compile(
            r'(?P<prefix>\w+:\d+:\w+:\d+:\d+:\d+:\d+(:[ACGTN]\+[ACGTN])? )'
            r'(?P<mate>[12])'
            r'(?P<suffix>:[YN]:\d*[02468]:([ACGTN]|\d)+)'
        )
    )
    # Illumina Casava <1.8
    # example: HWUSI-EAS100R:6:73:941:1973#0/1
    conventions.append(
        re.compile(
            r'(?P<prefix>[\w-]+:\d+:\d+:\d+:\d+#([ACGTN|\d])+/)'
            r'(?P<mate>[12])'
            r'(?P<suffix>)'
        )
    )

    # Scan for regexes and return mate information
    for regex in conventions:
        match = re.search(regex, seq_id)
        if match:
            mate_info = int(match.group('mate'))
            return (
                mate_info,
                ''.join([
                    match.group('prefix'),
                    match.group('suffix'),
                ])
            )

    # Return value if no regex matched
    return (0, seq_id)


def summarize_mate_info(
    counts: Dict[int, int]
) -> str:
    """Summarizes combined mate information for set of reads from counts for
    each mate.

    Args:
        counts: Dictionary of read counts for first and second mate, with keys
            `1` and `2` respectively.

    Returns:
        One of three values:
            "first_mate": All reads first mate
            "second_mate": All reads second mate
            "mixed_mates": First and second mate reads mixed

    Raises:
        ValueError: No mate counts.
    """
    # First mate
    if counts[1] and not counts[2]:
        return Outcomes.first_mate.value
    # Second mate
    if counts[2] and not counts[1]:
        return Outcomes.second_mate.value
    # First and second mates mixed
    if counts[1] and counts[2]:
        return Outcomes.mixed_mates.value
    # No reads
    raise ValueError("No records in file")
