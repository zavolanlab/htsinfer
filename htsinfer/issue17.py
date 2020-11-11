
sequences_DNA = [
    "CGTACTTGGTTCGATCAACCAGGCCGTAAAGTACGACGTCATCAAAAACGTTTAAACAAA",
    "GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC",
    "GCTGCTGGTCTTCATCCAACCTATGCTCGAACCATCGGTATTTCAGTTGATCATCGACGA",
    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"
    ]

sequences_RNA = [
    "CGUACUUGGUUCGAUCAACCAGGCCGUAAAGUACGACGUCAUCAAAAACGUUUAAACAAA",
    "GCUCAACGUGUUGCUCCUCGACCAGCCAAAGGUUCAUUACGGCCAGUUGUUCGUGGUACC",
    "GCUGCUGGUCUUCAUCCAACCUAUGCUCGAACCAUCGGUAUUUCAGUUGAUCAUCGACGA",
    "ACACGUCGUUAUAACAUGAAAGUACGUUCUGGACGCGGUUUUUCUUUGGAUGAAAUUCGU"
    ]

core_motifs_DNA = ["GCT", "TAT", "ACA", "CGC", "TGT", "GAG", "AAA", "GGG", "CCC", "TTT"]

core_motifs_RNA = ["GCU", "UAU", "ACA", "CGC", "UGU", "GUG", "AAA", "GGG", "CCC", "UUU"]

from typing import Dict
import math

def extend_motifs(input_sequences:list, core_motifs:list,
                  nucleic_acid:str):
    """ Extends motifs based on nucleotide representation in
        sequenced reads

        Args:
            input_sequences: list of sequences
            core_motifs: list of core motifs
            nucleic_acid: type of nucleic acid (dna, rna)

        Returns:
            extended_motifs: list of extended motifs

        Raises:
            TypeError: input_sequences is not a list
            TypeError: input_sequences is not a list
            ValueError: nucleic_acid is not in specified range

        """
"""
    if not isinstance(input_sequences, list):
        raise TypeError('input_sequences is not a list of sequences!')

    if not isinstance(core_motifs, list):
        raise TypeError('core_motifs is not a list of sequences!')

    if nucleic_acid.lower() not in ['rna', 'dna']:
        raise ValueError('nucleic_acid is not in specified range (dna/rna)!')


    base_dict = {
        'dna': 'ACTG',
        'rna': 'ACUG'
    }

    # determine valid characters depending on input sequences
    try:
        valid_bases = base_dict[nucleic_acid.lower()]
    except KeyError:
        raise ValueError(
            'nucleic_acid is not in specified range (dna/rna)!')
    return extended_motifs
"""

def find_motif_positions(read:str, coremotif:str):
    """
    Returns the start and end position of a core motif in a read.

    Args:
        read
        coremotif

    Returns:
        strartpositions: list of startpositions of core motif in read
        endpositions: list of endpositions of core motif in read
    """
    startpositions = [] #start position of core motif in sequence
    for i in range(len(read)-len(coremotif)+1): # loop over alignment
        match = True
        for j in range(len(coremotif)): # loop over characters
            if read[i+j] != t[j]:  # compare characters
                match = False   # mismatch
                break
        if match == True:   # allchars matched
            startpositions.append(i)
    print(startpositions)
    endpositions = []
    for k in range(len(startpositions)):
        endpositions.append(startpositions[k]+len(coremotif))
    print(endpositions)
    return startpositions, endpositions

s = "GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC"
t = "GCT"
start,end = find_motif_positions(s,t)


#find start and endpositions in all reads:
right_freq: Dict[str, int] = {'T': 0, 'C': 0, 'A': 0, 'G': 0}
left_freq: Dict[str, int] = {'T': 0, 'C': 0, 'A': 0, 'G': 0}
right_prob: Dict[str, int] = {}
left_prob: Dict[str, int] = {}
motif = "GCT"
for read in (input_sequences):
    start,end = find_motif_positions(read, motif)
    #iterating through all instances of the core motif in the read
    for pos in len(start):
        #left flanking position
        if start[pos]>0:
            char = read[start[pos]-1]
            freq = left_freq.get(char)
            dict_char = {char: freq+1}
            left_freq.update(dict_char)
        else:
            break
        #right flanking position
        if end[pos] < len(read):
            char = read[start[pos]-1]
            freq = right_freq.get(char)
            dict_char = {char: freq+1}
            right_freq.update(dict_char)
        else:
            break
    #left flanking position
    tot = sum(left_freq.values())
    nucleotides = ["T", "C", "A", "G"]
    entropy = 0
    for n in nucleotides:
        p = left_freq.get(n)/tot
        left_prob[n] = p
        entropy += p*(math.log2(p))

    #right flanking position
    tot = sum(left_freq.values())
    nucleotides = ["T", "C", "A", "G"]
    entropy = 0
    for n in nucleotides:
        p = left_freq.get(n)/tot
        left_prob[n] = p
        entropy += p*(math.log2(p))

    ##write functions to reduce redundancy!






