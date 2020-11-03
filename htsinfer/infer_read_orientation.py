"""Infer read orientation from sample data."""
import os, typing

from Bio import SeqIO


def infer():
    """Main function coordinating the execution of all other functions.
    Should be imported/called from main app and return results to it.
    """
    # implement me
    output_path = get_seq(filepath = filepath, organism = organism, temp_dir = temp_dir)


def get_transcripts(filepath: str, organism: typing.Union[str, int], temp_dir: str) -> str:
    """
    Selects transcript sequences of the desired organims,
    Example FASTA header: rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029
    organism short name: apisum, taxon id: 7029
    outputs it into a FASTA file

    Args:
        filepath(str): Path to the original FASTA file
        organism(Union[str, int]): organism short name(str, column 4 of FASTA header) or taxon id(int, column 5 of FASTA header)
        temp_dir(str): Path to the temporary directory

    Raises:
        OSError: if temp_dir does not exist

    Returns:
        output_path(str): Path to the output file saved in the temporary directory

    """
    output_path = os.path.join(temp_dir, organism + ".fasta")
    hits = []
    with open(filepath, "r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if (isinstance(organism, str) and record.description.split("|")[3] == organism) or (
                    isinstance(organism, int) and record.description.split("|")[4] == str(organism)):
                hits.append(record)
    SeqIO.write(hits, output_path, "fasta")
    return output_path
