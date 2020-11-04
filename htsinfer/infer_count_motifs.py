"""Infer motif frequency from sample data."""


from typing import Dict


def count_motifs(input_sequences: list, min_motif_length: int,
                 max_motif_length: int, nucleic_acid: str = 'dna',
                 non_nucleotide_characters: str = 'ignore_seqs') -> Dict:
    """Count frequency of motifs in one or more nucleotide sequences

        Args:
            input_sequences: list of sequences
            min_motif_length: minimal length of motif
            max_motif_length: maximal length of motif
            nucleic_acid: type of nucleic acid ('dna', 'rna')
            non_nucleotide_characters: handling of non-nucleotide chars
                ('ignore_chars', 'ignore_seqs', 'include')

        Returns:
            motif_freq: paired motif frequency {motif: frequency}

        Raises:
            TypeError: `input_sequences` is not a list
            ValueError: `nucleic_acid` is not in specified range
            ValueError: `min_motif_length` is not a positive integer
            ValueError: `max_motif_length` is not a positive integer
            ValueError: `min_motif_length` is longer than the shortest list
            ValueError: `min_motif_length` is longer than `max_motif_length`

    """

    if not isinstance(input_sequences, list):
        raise TypeError('input_sequences is not a list of sequences!')

    if not input_sequences:
        raise TypeError('input_sequences is not a list of sequences!')

    if min_motif_length <= 0:
        raise ValueError('min_motif_length is not a positive integer!')

    if max_motif_length <= 0:
        raise ValueError('max_motif_length is not a positive integer!')

    if min_motif_length > max_motif_length:
        raise ValueError('min_motif_length is longer than max_motif_length!')

    if min_motif_length > len(min(input_sequences)):
        raise ValueError('min_motif_length is longer than shortest sequence!')

    motif_freq: Dict[str, int] = {}  # Create empty dictionary (str -> int)

    for sequence_no in range(0, len(input_sequences)):
        seq = input_sequences[sequence_no]

        # determine valid characters depending on input sequences
        if nucleic_acid.lower() == 'dna':
            valid_bases = dict.fromkeys('ACTG')
        if nucleic_acid.lower() == 'rna':
            valid_bases = dict.fromkeys('ACUG')

        # non_nucleotide_characters == 'ignore_seqs' skips this sequence
        # and continues with the next sequence
        try:
            if not all(base in valid_bases for base in seq) \
                    and non_nucleotide_characters == 'ignore_seqs':
                print('sequence', sequence_no+1, 'ignored')
                continue
        except NameError:
            raise ValueError(
                'nucleic_acid is not in specified range (dna/rna)!')

        # collect all motifs and add to dict
        for motif_length in range(min_motif_length, max_motif_length + 1):
            for char in range(0, len(seq) - motif_length + 1):
                motif = seq[char: char + motif_length]
                if all(char in valid_bases for char in motif):
                    if motif in motif_freq:
                        motif_freq[motif] += 1
                    else:
                        motif_freq[motif] = 1
                elif non_nucleotide_characters == 'include':
                    if motif in motif_freq:
                        motif_freq[motif] += 1
                    else:
                        motif_freq[motif] = 1
                elif non_nucleotide_characters == 'ignore_chars':
                    print('motif not valid: ', motif)

    return motif_freq
