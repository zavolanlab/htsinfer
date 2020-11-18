'''Extends motifs based on nucleotide representation in sequenced reads'''

import math
import numpy as np
from scipy import stats

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
            if read[i+j] != coremotif[j]:  # compare characters
                match = False   # mismatch
                break
        if match == True:   # allchars matched
            startpositions.append(i)
    print(startpositions)
    endpositions = []
    for k in range(len(startpositions)):
        endpositions.append(startpositions[k]+len(coremotif)-1)
    print(endpositions)
    return startpositions, endpositions

#test
#s = "GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC"
#t = "GCT"
#start,end = find_motif_positions(s,t)

def compute_entropy(input_sequences: list, motif: str, position: str):
    """
    Returns the entropy at a flanking position of a core motif.

    Args:
        input_sequences
        motif
        position: left or right

    Returns:
        entropy: float
    """
    nucleotides = np.array(["T", "C", "A", "G"])
    numbers = np.array([0,1,2,3])
    freq = np.zeros(shape=4)
    flanking = np.empty(shape=len(input_sequences), dtype=str)
    for read_no, read in enumerate(input_sequences):
        start,end = find_motif_positions(read, motif)
        for pos in range(len(start)): #iterating through all instances of the core motif in the read
            if position == "left":
                if start[pos]>0:
                    char = read[start[pos]-1]
                    freq[[numbers[nucleotides == char]]] += 1
                    flanking[read_no] = char
                else:
                    break
            if position == "right":
                if end[pos] < (len(read)-1):
                    char = read[end[pos]+1]
                    freq[[numbers[nucleotides == char]]] += 1
                    flanking[read_no] = char
                else:
                    break
    flanking_mode = stats.mode(flanking)[0][0]
    print(flanking_mode)
    print(freq)
    tot = np.sum(freq)
    print(tot)
    entropy = 0
    for i in range(len(nucleotides)):
        p = freq[i]/tot + 0.00000000001
        entropy += p*(math.log2(p))
    print(entropy)
    return entropy, flanking_mode

#test
#motif = "AAG"
#position = "right"
#compute_entropy(sequences_DNA, motif, position)

def extend_motifs(input_sequences:list, core_motifs:list, nucleic_acid: str):
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
    if not isinstance(input_sequences, list):
        raise TypeError('input_sequences is not a list of sequences!')

    if not isinstance(core_motifs, list):
        raise TypeError('core_motifs is not a list of sequences!')

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
    for sequence_no, seq in enumerate(input_sequences):
        if not all(base in valid_bases for base in seq):
            raise ValueError('invalid bases for specified nucleic_acid in input_sequences')
    #code
    cutoff = 0.1
    extended_motifs = []
    for motif in core_motifs:
        position = "left"
        entropy, flanking_mode = compute_entropy(input_sequences, motif, position)
        extended_motif = motif
        while abs(entropy) < cutoff:
            extended_motif = flanking_mode + extended_motif
            print(extended_motif)
            entropy, flanking_mode = compute_entropy(input_sequences, extended_motif, position)
        position = "right"
        entropy, flanking_mode = compute_entropy(input_sequences, extended_motif, position)
        while abs(entropy) < cutoff:
            extended_motif = extended_motif + flanking_mode
            print(extended_motif)
            entropy, flanking_mode = compute_entropy(input_sequences, extended_motif, position)
        extended_motifs.append(extended_motif)
    print(extended_motifs)
    return extended_motifs

#test
#core_motifs = ["AAG", "AT", "CGT"]
#extend_motifs(sequences_DNA, core_motifs)







