#!/usr/bin/env python
# coding: utf-8

### Created on 04 November 2020 ###
###         Issue 21            ###
### Ji Hoon Han, Elifnaz Celik  ###

"""Infer potential barcode k-mers from sample data"""

from enum import Enum
from typing import List, Tuble


class Outcomes(Enum):
    invalid_nuc_requested = "invalid_nucleotide_requested"
    invalid_number_requested = "invalid_number_requested"

def group_kmers(input_kmers: List[str], kmer_size: float) -> List[str]:
    
    """Function that groups list of potential barcodes by kmer size
    
        Args:
            input_kmers (List[str]) : list of potential barcodes
            kmer_size(float) : K-mer length.
               
        Returns:
            selected_barcode (Tuple[str, float]): A list of tuples of selected barcode and the count {barcode, total count}
            
    
        Raises:
            TypeError:  'input_kmers' is not a list
            TypeError:  'kmer_size' is not a float
    
    """
    
    validate kmer character input parameters
    if kmer_input != ["A","C","T","G"]:
        return Outcomes.invalid_nuc_requested.value
    
    # Generate a list of selected barcodes of same kmer length
    selected_barcode = []
    
    # iteration for each barcode
    for barcode in kmer_input:
        # if the length barcode is same as the input
        if len(barcode) == kmer_size:
            # append the list
            selected_barcode.append(barcode)
    return selected_barcode, len(selected_barcode)


def similarity_kmers(barcode: List[str]) -> Tuples[str,str,int,float,str]: 
    
    """Function that returns position dependent frequency of list of grouped barcodes
    
        Args:
            selected_barcode (List[str]) : list of grouped barcodes
            
        Returns:
            barcodes (Tuple[str,str,int,float,str]) : A list of tuples of position-dependent frequency measured barcodes
                                                       i.e {a,b,c,d,e} where:
                                                       a = first barcode
                                                       b = second barcode
                                                       c = kmer_size
                                                       d = position-dependent frequency
                                                       e = position-dependent frequency score
        
        Raises:
            TypeError : 'selected_barcode' is not a list
    
    """
    
    loop_cnt = 0
    
    # Create a matrix for list of barcodes
    for barcode1 in barcode:
        loop_cnt2 = 0
        for barcode2 in barcode:
            if loop_cnt2 >= loop_cnt:
                score = 0
                # To remove redundancy, loop count is measured
                for i in range(0, len(barcode1)):
                    # Measuring the position-dependent frequency
                    if barcode1[i] == barcode2[i]:
                        score += 1
                    else:
                        score += 0
                # Score percentage by number of match (score) by the total length of kmer        
                score_percentage = float(score)/float(len(barcode1))
                
                # Categorising by score percentage but this can be modified by the user
                # Here categorised as "High score" for >= 0.7, "Low score" for <= 0.3,
                # and "Middle" for the rest
                if score_percentage >= 0.7:
                    score_value = "High score"
                elif score_percentage <= 0.3:
                    score_value = "Low score"
                else:
                    score_value = "Middle"
                
                # Logically, potential barcodes have higher abundance therefore
                # return list of barcodes with high position-dependent frequency score
                if score_value == "High score":
                    print(barcode1,barcode2,len(barcode1),score_percentage,score_value)
            loop_cnt2 += 1
        loop_cnt += 1
