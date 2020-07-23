"""Infer organism information for sequencing library."""

import os
import pandas as pd
import subprocess as sp


def kallisto(file_1: str, file_2: str = None):
    """Builds index file and run quantification algorithm.

    Args:
        file_1 (str) : File path to read/first mate library.
        file_2 (str) : File path to second mate library.
    """
    index = "kallisto index -i transcripts.idx --make-unique transcripts.fasta"
    quant_single = "kallisto quant -i transcripts.idx -o output -l 100 -s 300 --single " + file_1
    if file_2 is not None:
        quant_paired = "kallisto quant -i transcripts.idx -o output " + file_1 + " " + file_2

    try:
        with open(os.devnull, "w") as fnull:
            sp.call(index, shell=True)
            if file_2 is not None:
                sp.call(quant_paired, shell=True)
            else:
                sp.call(quant_single, shell=True)
    except sp.CalledProcessError:
        raise Exception("Error running kallisto")


def count_info():
    """Infers organisms count and return them as dictionary.

    Returns:
        A Dictionary with count percentage for every organism.
        example : ('scerevisiae', '4932'): 89.22159777956873
    """
    # Reading tsv file created by kallisto quant
    path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(path, "output\\abundance.tsv")
    df = pd.read_csv(file_path, sep='\t')

    # Dictionary to store organism info
    organism_count: Dict[(str, int), float] = {}
    dimension = df.shape
    rows = dimension[0]
    total_tpm = 0.0
    for i in range(rows):
        row_id = df['target_id'][i]
        contents = list(map(str, row_id.split("|")))
        organism_name = contents[3]
        organism_tax_id = contents[4]
        # Update organism count
        if (organism_name, organism_tax_id) in organism_count:
            organism_count[(organism_name, organism_tax_id)] += float(df['tpm'][i])
        else:
            organism_count[(organism_name, organism_tax_id)] = float(df['tpm'][i])
        total_tpm += float(df['tpm'][i])

    # Sorting as per organism with the highest counts
    sorted_organism_count = {k: v for k, v in sorted(organism_count.items(), key=lambda item: -1*item[1])}

    # Calculating Percentage
    for i in sorted_organism_count:
        sorted_organism_count[i] = (sorted_organism_count[i]/total_tpm)*100
    return sorted_organism_count
