## Task:
Given a set of sequences, calculate the number of occurrences of all motifs within a range of 
length that occur in the input sequences. Return the results as a dictionary.

## Function description:

1) for every sequence: determine all motifs/k-mers (k in range(a,b))  in the sequences with string slicing, save as list 	
2) make one big list combining all the smaller lists 
3) construct an output dictionary with motif sequences (key) and corresponding frequencies (value)
4) output sorted dictionary

## Function outline: 

def count_motifs(input_sequences: list, min_motif_length: int, max_motif_length: int) -> dict
""" Function that calculates the occurrence of all possible motifs in one or multiple sequences and 
    returns a dictionary with all motifs within the specified length and their occurrence

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