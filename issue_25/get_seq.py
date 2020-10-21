# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 11:25:08 2020

@author: Bátora Dániel
"""

import pickle as pkl
import os 
from Bio import SeqIO, Seq
from typing import Union



def get_seq(filepath:str,organism:Union[str, int],temp_dir:str) -> str:
    
    """
    Selects DNA sequence of the desired organims,
    outputs it into a FASTA file 
    
    Args: 
        filepath(str): Path to the original FASTA file 
        organism(Union[str, int]): organism short name(str) or taxon id(int)
        temp_dir(str): Path to the temporary directory
        
    Returns: 
        output_path(str): Path to the output file saved in the temporary directory
    
    """
    
    
    
    if not os.path.isdir(temp_dir): 
        os.mkdir(temp_dir)
    output_path = os.path.join(temp_dir, organism + ".fasta")
    reads = []
    with open(filepath,  "r", encoding = "utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            reads.append(record)
    
    hits = []
    for record in reads:
        if type(organism) == str: 
            if record.description.split("|")[3] == organism:
                hits.append(record)
        elif type(organism) == int: 
            if record.description.split("|")[4] == str(organism):
                hits.append(hits)
    SeqIO.write(hits, output_path, "fasta")
    
    return output_path
        


hits = get_seq(os.path.join(os.getcwd(), "transcripts.fasta"),  "avaga", "C:/github/htsinfer/data/temp_dir")



