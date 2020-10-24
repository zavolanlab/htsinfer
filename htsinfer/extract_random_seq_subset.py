"""Extract a random subset of the sample data."""

import shutil
import random
import gzip
import logging
from typing import List, Tuple, IO
from enum import Enum
from Bio.SeqIO.QualityIO import FastqGeneralIterator  # type: ignore


class Outcomes(Enum):
    """Enumerator for processing outcomes."""
    files_verified = "files_verified"
    no_input_provided = "no_input_provided"
    no_output_provided = "no_output_provided"
    no_input_output_provided = "no_input_output_provided"
    mismatched_in_out_files = "mismatched_in_out_files"
    same_in_out_files = "same_in_out_files"
    extraction_done = "extraction_done"
    file_error = "file_error"
    no_input_sequences = "no_input_sequences"
    invalid_number_requested = "invalid_number_requested"
    unexpected_return_value = "unexpected_return_value"


LOGGER = logging.getLogger(__name__)


def extract(
    infile: List[str],
    outfile: List[str],
    proportion: float = 1.0,
        max_records: int = 10000) -> str:

    """Extract a subset of reads from one or two FASTQ files.

    Args:
        infile (List[str]) : List of file paths to read libraries
        outfile (List[str]) : List of file paths to extracted reads
        proportion (float) : Proportion of total number of
            sequences to extract, if not exceeding max_records
            and not exceeding the number of input sequences.
        max_records: Maximum number of reads to extract, if not
            exceeding number of input sequences.
            Set to `0` to extract all records.

    Returns:
        outcome of processing (str).

    Examples:
        >>> extract(
        ...     infile=["../tests/test_files/first_mate.fastq"],
        ...     outfile=["../tests/test_files/first_mate_subset.fastq"],
        ...     max_records=10)

        >>> extract(
        ...     infile=["../tests/test_files/first_mate.fastq",
                        "../tests/test_files/second_mate.fastq"],
        ...     outfile=["../tests/test_files/first_mate_subset.fastq",
                        "../tests/test_files/second_mate_subset.fastq"],
        ...     max_records=10)
    """

    # validate numerical input parameters
    if proportion <= 0.0 or max_records < 0:
        return Outcomes.invalid_number_requested.value

    # check that we got valid files as parameters
    check_files = verify_files(infile, outfile)
    if check_files != Outcomes.files_verified.value:
        LOGGER.debug(f"Invalid file arguments error: {check_files}")
        return check_files

    # Process files
    for (inputfile, outputfile) in zip(infile, outfile):
        LOGGER.debug(f"Processing file: {inputfile}")
        outcome = process_fastq_file(
            infile=inputfile,
            outfile=outputfile,
            proportion=proportion,
            max_records=max_records)
        if not outcome == Outcomes.extraction_done.value:
            LOGGER.debug(
                f"Error trying to extract from {inputfile} to {outputfile}")
            return outcome

    LOGGER.debug("Done extracting sequences...")
    return Outcomes.extraction_done.value


def process_fastq_file(
    infile: str,
    outfile: str,
    proportion: float = 1.0,
        max_records: int = 10000) -> str:

    """Extracts a  subset of sequences from a single FASTQ file.

    Args:
        infile (str): FASTQ file to process.
        outfile (str): FASTQ file of extracted sequences.
        proportion (float): Proportion of sequences to extract
            if not exceeding max_records.
        max_records (int): Maximum number of records to extract.
        Set to `0` to extract all records.

    Returns:
        status (Outcomes)
    """

    # validate input file
    (input_lines, status) = count_lines(infile)

    if status != "files_verified":
        return status

    # we got input
    # deal with the base case first:
    # if max_records set to 0, all sequences should be extracted
    # so we just make a copy of the file
    if max_records == 0:
        LOGGER.debug("Copying full file(s)...")
        shutil.copy2(infile, outfile)
        return Outcomes.extraction_done.value

    # we have to extract a subset of sequences
    # we know how many lines we have in the input file
    # number of sequences will be lines/4 lines (per seq)
    lines_per_record = 4  # given by the FASTQ format

    seqs_ids = choose_seq_indices(
            int(input_lines/lines_per_record),
            proportion,
            max_records)

    # we read the input file
    # when we reach an index from seqs_ids, we save the record

    LOGGER.debug("Opening file...")

    # open input file
    try:
        if infile.endswith(".gz"):
            fin = gzip.open(infile, 'rt')
        else:
            fin = open(infile, 'rt')
    except OSError:
        LOGGER.debug(
            f"Invalid input file '{infile}'. Error: Could not open file")
        return Outcomes.file_error.value

    # open output file
    try:
        with open(outfile, 'wt') as fout:
            # save selected sequences
            outcome = write_out_seqs(fin, fout, seqs_ids)
        # done
        # close input file
        fin.close()
        return outcome

    except OSError:
        LOGGER.debug(
            f"Invalid input file '{outfile}'. Error: Could not open file")
        return Outcomes.file_error.value


def write_out_seqs(
    fin: IO,
    fout: IO,
        seqs_ids: List[int]) -> str:

    """Extracts a  subset of sequences from a file handle
        and writes them to another file handle.

    Args:
        fin (file handle): FASTQ file object to process.
        fout (file handle): FASTQ file object of extracted sequences.
        proportion: (float) Proportion of sequences to extract
            if not exceeding max_records and number of input sequences.
        seqs_ids (List(int)): indices of sequences from input
            file object to write to output file object
    Returns:
        Status (Outcomes)
    """

    # keep track of the index of the read in the input file
    # keep a the pointer in seqs_ids
    ind = 0
    records = -1

    # loop through the read records
    try:
        for (seqid, seq, qual) in FastqGeneralIterator(source=fin):
            # stop if we have extracted all the records we wanted
            if ind >= len(seqs_ids):
                break
            # update the index of the current record
            records += 1
            # check if this is a record we should write out
            if records == seqs_ids[ind]:
                # write out the record
                fout.write("@" + seqid + "\n")
                fout.write(seq + "\n")
                fout.write("+\n" + qual + "\n")

                # move the pointer in the list seqs_to_get
                ind += 1
    except ValueError:
        return Outcomes.no_input_sequences.value

    return Outcomes.extraction_done.value


def count_lines(infile: str) -> Tuple[int, str]:
    """Count the lines in a (gzipped) FASTQ file.

    Args:
        infile (str): FASTQ file to process.

    Returns:
        nr_lines (int)
        status (Outcomes)
    """

    #   count the lines in a (gzipped) FASTQ file
    lines = 0

    try:
        if infile.endswith(".gz"):
            with gzip.open(infile, 'r') as fin:
                for _ in fin:
                    lines += 1
        else:
            with open(infile, 'r') as fin:
                for _ in fin:
                    lines += 1
    except OSError:
        LOGGER.debug(
            f"Invalid input file '{infile}'. Error: Could not open file")
        return (lines, Outcomes.file_error.value)

    if lines:
        return (lines, Outcomes.files_verified.value)

    return (lines, Outcomes.no_input_sequences.value)


def choose_seq_indices(
    input_seqs: int,
    proportion: float = 1.0,
        max_records: int = 10000) -> List[int]:

    """Generate a list seq indices.
    Args:
        input_seq (int): number of input_sequences
        proportion (float): proportion of seqs to extract
            if not larger than max_records and not larger
            than nr of input sequences
        max_records (int): max_nr of records to extract if
            not larger than the number of input sequences

    Returns:
        list of indices (list(int))
    """

    # default number of sequences to get: entire set
    nr_seqs_to_get = input_seqs

    # should we extract a specified proportion?
    if proportion < 1.0:
        nr_seqs_to_get = int(proportion * input_seqs)
    elif max_records < nr_seqs_to_get:
        # we want a specific number
        # if it's not larger than the number of sequences we have
        nr_seqs_to_get = max_records

    # get a random subset of sequence indices
    # sorted from lowest to highest
    seqs_ids = list(range(input_seqs))
    random.shuffle(seqs_ids)

    # get the indices of sequences to extract, sorted
    seqs_ids = seqs_ids[:nr_seqs_to_get]

    return sorted(seqs_ids)


def verify_files(
    infile: List[str],
        outfile: List[str]) -> str:

    """Error checking on input/output files

    Args:
        infile (List[str]) : List of file paths to read libraries
        outfile (List[str]) : List of file paths to extracted reads

    Returns:
        status (Outcomes).
    """

    ret_value = ""
    if not len(infile) > 0:
        LOGGER.debug("No input files provided...")
        ret_value = Outcomes.no_input_provided.value

    if not len(outfile) > 0:
        LOGGER.debug("No output files provided..")
        if ret_value:
            ret_value = Outcomes.no_input_output_provided.value
        else:
            ret_value = Outcomes.no_output_provided.value

    if ret_value:
        return ret_value

    if not len(infile) == len(outfile):
        LOGGER.debug(
            "Numbers of input and output files don't match..")
        return Outcomes.mismatched_in_out_files.value

    for (inputfile, outputfile) in zip(infile, outfile):
        if inputfile == outputfile:
            LOGGER.debug("Output file cannot be the same as the input file...")
            return Outcomes.same_in_out_files.value

    return Outcomes.files_verified.value
