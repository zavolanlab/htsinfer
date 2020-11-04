"""Infer potential barcode k-mers"""

import pandas as pd
from enum import Enum
from typing import List, Dict


class Outcomes(Enum):
    invalid_number_requested = "invalid_number_requested"

def kmer_similarity(kmer_input: list, kmer_size: int >= 5 or <= 20) -> Dict:
    
    """Function that returns set of kmers
    
        Args:
            kmer_input (list) : list of k-mers
            kmer_size (integer) : length of k-mers, minimum 5 and maximum 20
                    
        Returns:
            set of kmers
            
        Raises:
            TypeError:  'input_kmers' is not a list
            TypeError:  'kmer_size' is not a integer
            ValueError: 'kmer_size' value is invalid
   
    """

    # validate numerical input parameters
    if kmer_size <= 0:
        return Outcomes.invalid_number_requested.value
       
    set_of_kmers = Dict[str] = {}
    
    for kmer_no in len(kmer_input):
        kmer = kmer_input[kmer_no]
        
        kmer_set = kmer[]
        

        
    return set_of_kmers


"""
Example

Potential_barcode=["ACT","TTTCT","TCTACT","TCTT","GCGTTTTAT"]
test = kmer_similarity(kmer_input=Potential_barcode, kmer_size=3)
print(test)
"""
