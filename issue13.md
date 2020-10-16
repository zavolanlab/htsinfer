## Task:
Given a set of sequences, calculate the number of occurrences of all motifs within a range of 
length that occur in the input sequences. Return the results as a dictionary.

## Function description:

input parameters:
* sequences = list of sequences
* a,b: specifying the range of lengths of motifs 
output: dictionary with paired data: { "motif_seq" : frequency }

1) create all possible motifs with lengths in range [a,b], save as list motifs
idea for code:
nucleotides = ["A", "T", "G", "C"]
for i in range(a,b):
	motifs = itertools.product(nucleotides, repeat=i)
already add these in a dict with all freq counters = 0
	
2) scan all sequences for all motifs and record frequencies
e.g. string.count(substring) --> sequence.count(motif)
BUT: "the needle in a haystack problem" --> runtime...
fastest method might be regex: 
for motif in motifs:
    if motif in sequence:
        dict[motif] += 1
-> see https://stackoverflow.com/questions/4901523/whats-a-faster-operation-re-match-search-or-str-find

3) construct and output dictionary with motif sequences and corresponding frequencies

## Function outline: 

def count_motifs(input_sequences: list, min_motif_length: int, max_motif_length: int) -> dict
""" Function that calculates the occurrence of all possible motifs in one or multiple sequence and 
    returns a dict with all motifs within the specified length and their occurrence

Args: 
    input_sequences (list): list of sequences
    min_motif_length (int): minimal length of motif
    max_motif_length (int): maximal length of motif
   
Returns:
    motif frequency (dict) with paired data {"motif_seq": frequency }
    
Raises:
    ValueError: input_sequences contains non-standard (ACTG) characters
    TypeError: input_sequences is not a list
"""