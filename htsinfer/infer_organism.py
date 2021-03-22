"""Infer organism information for sequencing library."""

import logging
import os
import subprocess as sp
from typing import (Dict, Tuple)
import zipfile as zp

import pandas as pd  # type: ignore

from htsinfer.utils import minmatch_factor_validator

LOGGER = logging.getLogger(__name__)


def infer(
    transcript_fasta: str,
    file_1: str,
    file_2: str = None,
    min_match: float = 10,
    factor: float = 2,
) -> str:
    """Builds index file and run quantification algorithm.

    Args:
        file_1: File path to first mate library.
        file_2: File path to second mate library.
        min_match: Minimum percentage that given organism needs to have
            to be considered as the resulting organism. If no organism is
            found more frequently than the specified number (in percent),
            null is returned in the JSON result to indicate that no organism
            could be confidently identified.
        factor: The minimum frequency ratio between the first and second most
            frequent organism in order for an organism sequence to be returned
            in the JSON result. If frequency ratio is less than the specified
            value, null is returned in the JSON result to indicate that no
            single organism could be confidently identified.

    Returns:
        Organism name.
    """
    with zp.ZipFile(transcript_fasta, "r") as zip_ref:
        zip_ref.extractall()

    # Runs Kallisto index
    exit_code = kallisto_index()
    if exit_code:
        return "NA"
    # Runs Kallisto quant
    exit_code = kallisto_quant(file_1=file_1, file_2=file_2)
    if exit_code:
        return "NA"

    LOGGER.debug("Processing organism count info")
    organism_tpm_count = process_count_info()
    organism_df = convert_dic_to_df(organism_tpm_count)
    # Checking confidence score
    if minmatch_factor_validator(organism_df, 0, min_match, factor):
        result = organism_df.iloc[0]['Organism']
    else:
        result = "NA"
    return result


def kallisto_index() -> int:
    """Builds an index from a FASTA formatted file of target sequences.

    Returns:
        Exit status of subprocess.
    """
    LOGGER.debug("Running Kallisto index")
    index_cmd = "kallisto index -i transcripts.idx --make-unique" + \
        " transcripts.fasta"
    result = sp.run(index_cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        LOGGER.debug(result.stderr)
    else:
        LOGGER.error(result.stderr)
    return result.returncode


def kallisto_quant(
    file_1: str,
    file_2: str = None,
) -> int:
    """Runs the quantification algorithm.

    Args:
        file_1: File path to first mate library.
        file_2: File path to second mate library.

    Returns:
        Exit status of subprocess.
    """
    LOGGER.debug("Running Kallisto quant")
    if file_2 is not None:
        quant_cmd = "kallisto quant -i transcripts.idx -o output " + \
            file_1 + " " + file_2
    else:
        quant_cmd = "kallisto quant -i transcripts.idx -o output" + \
            " -l 100 -s 300 --single " + file_1
    result = sp.run(quant_cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        LOGGER.debug(result.stderr)
    else:
        LOGGER.error(result.stderr)
    return result.returncode


def process_count_info() -> Dict[Tuple[str, int], float]:
    """Process organisms count info and return them as dictionry.

    Returns:
        A dictionary with count percentage for all organisms.
    """
    # Reading tsv file created by kallisto quant
    path = os.getcwd()
    path = os.path.join(path, "output", "abundance.tsv")
    abundance_df = pd.read_csv(path, sep='\t')

    # Dictionary to store organism info
    organism_tpm_count: Dict[Tuple[str, int], float] = {}
    dimension = abundance_df.shape
    rows = dimension[0]
    total_tpm: float = 0.0
    for i in range(rows):
        row_id = abundance_df['target_id'][i]
        contents = list(map(str, row_id.split("|")))
        organism_name = contents[3]
        organism_tax_id = int(contents[4])
        # Update organism tpm count
        if (organism_name, organism_tax_id) in organism_tpm_count:
            organism_tpm_count[(organism_name, organism_tax_id)] += \
                float(abundance_df['tpm'][i])
        else:
            organism_tpm_count[(organism_name, organism_tax_id)] = \
                float(abundance_df['tpm'][i])
        total_tpm += float(abundance_df['tpm'][i])

    # Calculating Percentage
    if total_tpm != 0:
        for (_name, _id), _ in organism_tpm_count.items():
            organism_tpm_count[(_name, _id)] = round(
                (organism_tpm_count[(_name, _id)]/total_tpm)*100, 2
                )
    return organism_tpm_count


def convert_dic_to_df(
    organism_tpm_count: Dict[Tuple[str, int], float]
) -> pd.DataFrame:
    """Converting dictionary into dataframe and writing json file.

    Args:
        organism_tpm_count: Dictionary of organism name, taxon id and
        it's count percentage.

    Retunrns:
        Dataframe of count info percentage for all organisms.
    """
    organism_df = pd.DataFrame(organism_tpm_count.items())
    organism_df[['Organism', 'Taxon ID']] = pd.DataFrame(
        organism_df[0].tolist()
        )
    organism_df = organism_df.sort_values(
        by=1, ascending=False
        ).reset_index(drop=True).drop([0], axis=1)
    organism_df = organism_df.rename(columns={1: 'Match %'})
    LOGGER.debug("Creating organism_count_info.json")
    organism_df.to_json(
        'organism_count_info.json', orient='split', index=False, indent=True
        )
    return organism_df
