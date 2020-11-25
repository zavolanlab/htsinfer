""" Extends motifs based on nucleotide representation in sequenced reads """

import math
import numpy as np  # type: ignore
from scipy import stats  # type: ignore

sequences_DNA = [
    "CGTACTTGGTTCGATCAACCAGGCCGTAAAGTACGACGTCATCAAAAACGTTTAAACAAA",
    "GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC",
    "GCTGCTGGTCTTCATCCAACCTATGCTCGAACCATCGGTATTTCAGTTGATCATCGACGA",
    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"
    ]


def find_motif_positions(sequence: str, motif: str):
    """
    Returns the start and end position(s) of a core motif in a sequence

    Args:
        sequence: string of nucleotides
        motif: string of nucleotides representing a motif

    Returns:
        startpositions: list of start position(s) of core motif in read
        endpositions: list of end position(s) of core motif in read
    """
    startpositions = []
    for i in range(len(sequence) - len(motif) + 1):  # loop over read
        match = True
        for j in range(len(motif)):  # loop over characters
            if sequence[i + j] != motif[j]:  # compare characters
                match = False   # mismatch
                break
        if match:   # all chars matched
            startpositions.append(i)
    endpositions = []
    for k in range(len(startpositions)):
        endpositions.append(startpositions[k] + len(motif) - 1)
    return startpositions, endpositions


def compute_entropy(input_sequences: list, motif: str, position: str):
    """
    Returns the entropy at a flanking position (left or right)
    of a motif in a list of sequences

    Args:
        input_sequences: list of sequences
        motif: string of nucleotides representing a motif
        position: left or right flanking position ('left', 'right')

    Returns:
        entropy: Shannon entropy (float) at flanking position
    """
    nucleotides = np.array(["T", "C", "A", "G"])
    numbers = np.array([0, 1, 2, 3])
    freq = np.zeros(shape=4)
    flanking = np.empty(shape=len(input_sequences), dtype=str)
    for seq_no, sequence in enumerate(input_sequences):
        start, end = find_motif_positions(sequence, motif)
        for pos in range(len(start)):
            if position == "left":
                if start[pos] > 0:
                    char = sequence[start[pos]-1]
                    freq[[numbers[nucleotides == char]]] += 1
                    flanking[seq_no] = char
                else:
                    break
            if position == "right":
                if end[pos] < (len(sequence)-1):
                    char = sequence[end[pos]+1]
                    freq[[numbers[nucleotides == char]]] += 1
                    flanking[seq_no] = char
                else:
                    break
    flanking_mode = stats.mode(flanking)[0][0]
    tot = np.sum(freq)
    entropy = 0
    for i in range(len(nucleotides)):
        prob = freq[i]/tot + 0.00000000001
        entropy += prob*(math.log2(prob))
    return entropy, flanking_mode


def extend_motifs(sequences: list, motifs: list, nucleic_acid: str):
    """ Extends motifs based on nucleotide representation in
        sequenced reads

        Args:
            sequences: list of sequences
            motifs: list of motifs
            nucleic_acid: type of nucleic acid ('dna', 'rna')

        Returns:
            extended_motifs: list of extended motifs

        Raises:
            TypeError: sequences is not a list
            TypeError: motifs is not a list
            ValueError: nucleic_acid is not in specified range
            ValueError: longest motif in motifs longer than shortest sequence
            ValueError: sequences or motifs is not a list
            ValueError: empty string in sequences or motifs

        """
    if not isinstance(sequences, list):
        raise TypeError('sequences is not a list!')

    if not isinstance(motifs, list):
        raise TypeError('motifs is not a list!')

    if len(max(motifs)) > len(min(sequences)):
        raise ValueError('longest motif in motifs '
                         'longer than shortest sequence!')

    if len(sequences) == 0 or len(motifs) == 0:
        raise ValueError('empty list')

    if len(min(sequences)) == 0 or len(min(motifs)) == 0:
        raise ValueError('empty string')

    base_dict = {
        'dna': 'ACTG',
        'rna': 'ACUG'
    }
    # determine valid characters depending on sequences
    try:
        valid_bases = base_dict[nucleic_acid.lower()]
    except KeyError:
        raise ValueError(
            'nucleic_acid is not in specified range (dna,rna)!')
    for seq_no, sequence in enumerate(sequences):
        if not all(base in valid_bases for base in sequence):
            raise ValueError('invalid bases for specified '
                             'nucleic_acid in sequences')

    cutoff = 0.1
    extended_motifs = []
    for motif in motifs:
        position = "left"
        entropy, flanking_mode = compute_entropy(sequences, motif, position)
        extended_motif = motif
        while abs(entropy) < cutoff:
            extended_motif = flanking_mode + extended_motif
            entropy, flanking_mode = compute_entropy(sequences,
                                                     extended_motif, position)
        position = "right"
        entropy, flanking_mode = compute_entropy(sequences,
                                                 extended_motif, position)
        while abs(entropy) < cutoff:
            extended_motif = extended_motif + flanking_mode
            entropy, flanking_mode = compute_entropy(sequences,
                                                     extended_motif, position)
        extended_motifs.append(extended_motif)
    return extended_motifs
