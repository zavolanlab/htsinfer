# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 11:25:08 2020

@author: Bátora Dániel
"""


import os 
from Bio import SeqIO
folder = "C:\github\htsinfer\data"
os.chdir(folder)
reads = []
with open("transcripts.fasta", "rU", encoding = "utf-8") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        reads.append(record)



def get_seq(file, organism):
    """
    Returns all sequences of an organism into a dict 

    Parameters
    ----------
    file : list
        List of SeqRecord objects, each containing 
        information about the sequence and the sequence .
    organism : str, int
        Organism to get

    Returns
    -------
    data : dict
        All instances of the organism in a dict with key being the organism
        and the value is a string containing the sequence.

    """
    assert type(file) == list
    assert type(organism) == str or type(organism) == int
    data = {}
    counter = 0 
    for record in file:
        if type(organism) == str: 
            if record.description.split("|")[3] == organism:
                data.update({organism + "_" + str(counter) : str(record.seq)}) 
                counter += 1  
        elif type(organism) == int: 
            if record.description.split("|")[4] == str(organism):
                data.update({str(organism) + "_" + str(counter) : str(record.seq)}) 
                counter += 1  
            
    return data
        
   