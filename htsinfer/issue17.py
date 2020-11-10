
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

def find_motif_positions(read, coremotif):
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
