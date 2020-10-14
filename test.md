Function description:

input parameters:
* sequences = list of sequences
* a,b: specifying the range of lengths of motifs 
output: dictionary with paired data: { "motif_seq" : frequency }

1) create all possible motifs with lengths in range [a,b], save as list motifs
idea for code:
nucleotides = ["A", "T", "G", "C"]
for i in range(a,b):
	motifs = itertools.product(nucleotides, repeat=i)
2) scan all sequences for all motifs and record frequencies
e.g. string.count(substring) --> sequence.count(motif)
BUT: "the needle in a haystack problem" --> runtime...
3) construct and output dictionary with motif sequences and corresponding frequencies
