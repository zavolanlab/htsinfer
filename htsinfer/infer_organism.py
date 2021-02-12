"""Infer organism information for sequencing library."""

import os
import pandas as pd
import subprocess as sp
import operator
import logging
from typing import (Dict, Tuple)

logger = logging.getLogger(__name__)


def kallisto(file_1: str, file_2: str = None) -> pd.DataFrame:
    """Builds index file and run quantification algorithm.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.

    Returns:
        A dataframe of count percentage information for top five organisms.
    """
    index = "kallisto index -i transcripts.idx --make-unique transcripts.fasta"
    quant_single = "kallisto quant -i transcripts.idx -o output -l 100 -s 300 --single " + file_1
    if file_2 is not None:
        quant_paired = "kallisto quant -i transcripts.idx -o output " + file_1 + " " + file_2

    try:
        with open(os.devnull, "w") as fnull:
            logger.debug("Running Kallisto index")
            sp.call(index, shell=True)
            logger.debug("Running Kallisto quant")
            if file_2 is not None:
                sp.call(quant_paired, shell=True)
            else:
                sp.call(quant_single, shell=True)
    except sp.CalledProcessError:
        logger.error("Invalid input file")
        raise Exception("Error : running kallisto/transcripts.fasta not located")

    logger.debug("Processing organism count info")
    result = process_count_info()
    # Converting result into dataframe
    oragnism_df = pd.DataFrame(result.items())
    oragnism_df.columns = ['Organism, Taxon ID', 'Match %']
    oragnism_df.index = oragnism_df.index + 1
    return oragnism_df


def process_count_info() -> Dict[Tuple[str, str], float]:
    """Infers organisms count and return them as dictionary.

    Returns:
        A Dictionary with count percentage for organisms.
        example : ('scerevisiae', '4932'): 89.22
    """
    # Reading tsv file created by kallisto quant
    path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(path, "output", "abundance.tsv")
    df = pd.read_csv(file_path, sep='\t')

    # Dictionary to store organism info
    organism_tpm_count: Dict[(str, int), float] = {}
    dimension = df.shape
    rows = dimension[0]
    total_tpm = 0.0
    for i in range(rows):
        row_id = df['target_id'][i]
        contents = list(map(str, row_id.split("|")))
        organism_name = contents[3]
        organism_tax_id = contents[4]
        # Update organism tpm count
        if (organism_name, organism_tax_id) in organism_tpm_count:
            organism_tpm_count[(organism_name, organism_tax_id)] += float(df['tpm'][i])
        else:
            organism_tpm_count[(organism_name, organism_tax_id)] = float(df['tpm'][i])
        total_tpm += float(df['tpm'][i])

    # Calculating Percentage
    for i in organism_tpm_count:
        organism_tpm_count[i] = round((organism_tpm_count[i]/total_tpm)*100, 2)

    # Sorting as per organism with the highest counts
    # sorted_organism_count = {k: v for k, v in sorted(organism_count.items(), key=lambda item: -1*item[1])}

    # Returning top five organisms as per TPM counts
    top_five = dict(sorted(organism_tpm_count.items(), key=operator.itemgetter(1), reverse=True)[:5])
    return top_five
