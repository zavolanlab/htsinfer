"""Infer organism information for sequencing library."""

import os
import pandas as pd
import subprocess as sp
import logging
import zipfile as zp
from typing import Dict

logger = logging.getLogger(__name__)


def kallisto(
    file_1: str,
    file_2: str = None,
    min_match: float = 10,
    factor: float = 2
) -> str:
    """Builds index file and run quantification algorithm.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.
        min_match (float) : minimum match percentage that top organism
        needs to have.
        factor (float) : factor by which first organism is greater than
        the second.

    Returns:
        Organism name if it satisfies confidence score else NA.
    """
    path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(path, "transcripts.fasta.zip")
    with zp.ZipFile(file_path, "r") as zip_ref:
        zip_ref.extractall()
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
    organism_df = process_count_info()
    # Checking confidence score
    if(confidence(organism_df, min_match, factor)):
        result = organism_df.iloc[0]['Organism']
    else:
        result = "NA"
    return result


def process_count_info() -> pd.DataFrame:
    """Infers organisms count info and return them as dataframe.

    Returns:
        A dataframe with count percentage for all organisms.
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
    if total_tpm != 0:
        for i in organism_tpm_count:
            organism_tpm_count[i] = round((organism_tpm_count[i]/total_tpm)*100, 2)

    # Sorting as per organism with the highest counts
    # sorted_organism_count = {k: v for k, v in sorted(organism_count.items(), key= lambda item: -1*item[1])}

    # Converting dictionary into dataframe
    organism_df = pd.DataFrame(organism_tpm_count.items())
    organism_df[['Organism', 'Taxon ID']] = pd.DataFrame(organism_df[0].tolist())
    organism_df = organism_df.sort_values(by=1,ascending=False).reset_index(drop=True).drop([0], axis=1)
    organism_df = organism_df.rename(columns={1: 'Match %'})
    organism_df.to_csv('organism_count_info.csv')
    return organism_df


def confidence(
    organism_df: pd.DataFrame,
    min_match: float = 10,
    factor: float = 2
) -> bool:
    """Checks confidence score

    Args:
        organism_df (dataframe) : count info percentage for all organisms,
        min_match (float) : minimum match percentage that top organism needs
        to have.
        factor (float) : factor by which first organism is greater than the
        second.

    Returns:
        bool : whether it satisfies confidence score or not.
    """
    if organism_df.iloc[0]['Match %'] < min_match:
        return False
    if organism_df.iloc[1]['Match %'] != 0:
        ratio = organism_df.iloc[0]['Match %']/organism_df.iloc[1]['Match %']
        if ratio < factor:
            return False
    return True
