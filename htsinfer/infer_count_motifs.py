"""Infer motif frequency from sample data."""


from typing import Dict


def count_motifs(input_sequences: list, min_motif_length: int,
                 max_motif_length: int, nucleic_acid: str) -> Dict:
    """Count frequency of motifs in one or more nucleotide sequences

        Args:
            input_sequences: list of sequences
            min_motif_length: minimal length of motif
            max_motif_length: maximal length of motif
            nucleic_acid: type of nucleic acid ('dna', 'rna')

        Returns:
            sorted_dict: paired motif frequency {motif: frequency}

        Raises:
            TypeError: `input_sequences` is not a list
            TypeError: `min_motif_length` is not an integer
            TypeError: `max_motif_length` is not an integer
            ValueError: `nucleic_acid` is not in specified range
            ValueError: `min_motif_length` is not a positive integer
            ValueError: `max_motif_length` is not a positive integer
            ValueError: `max_motif_length` is longer than the shortest list
            ValueError: `min_motif_length` is longer than the shortest list
            ValueError: `min_motif_length` is longer than `max_motif_length`

    """

    if not isinstance(input_sequences, list):
        raise TypeError('input_sequences is not a list of sequences!')

    if not input_sequences:
        raise TypeError('input_sequences is not a list of sequences!')

    if not isinstance(min_motif_length, int):
        raise TypeError('min_motif_length is not an integer!')

    if not isinstance(max_motif_length, int):
        raise TypeError('max_motif_length is not an integer!')

    if min_motif_length <= 0:
        raise ValueError('min_motif_length is not a positive integer!')

    if max_motif_length <= 0:
        raise ValueError('max_motif_length is not a positive integer!')

    if min_motif_length > max_motif_length:
        raise ValueError('min_motif_length is longer than max_motif_length!')

    if max_motif_length > len(min(input_sequences)):
        raise ValueError('max_motif_length is longer than shortest sequence!')

    if min_motif_length > len(min(input_sequences)):
        raise ValueError('min_motif_length is longer than shortest sequence!')

    if nucleic_acid.lower() not in ['rna', 'dna']:
        raise ValueError('nucleic_acid is not in specified range (dna/rna)!')

    motif_freq: Dict[str, int] = {}  # Create empty dictionary (str -> int)

    for sequence_no in range(0, len(input_sequences)):
        seq = input_sequences[sequence_no]
        user_choice = "0"

        # determine valid characters depending on input sequences
        if nucleic_acid.lower() == 'dna':
            valid_bases = dict.fromkeys('ACTG')
        if nucleic_acid.lower() == 'rna':
            valid_bases = dict.fromkeys('ACUG')

        # check if all characters in current sequence are valid.
        # if not, ask user what to do with invalid characters
        if not all(base in valid_bases for base in seq):
            print('sequence', sequence_no+1, 'contains unknown characters')
            print('you have the following choices:')
            print('(1) ignore the non-nucleotide-characters'
                  'and exclude them from the motifs')
            print('(2) include the non-nucleotide characters')
            print('(3) ignore the sequence with the non-nucleotide characters')
            user_choice = input('please specify how you want'
                                'to proceed using the numbers above: ')

        # user choice 3 skips this sequence
        # and continues with the next sequence
        if user_choice == '3':
            print('sequence', sequence_no+1, 'ignored')
            continue

        # collect all motifs and add to dict
        for motif_length in range(min_motif_length, max_motif_length + 1):
            for char in range(0, len(seq) - motif_length + 1):
                motif = seq[char: char + motif_length]
                if all(char in valid_bases for char in motif):
                    if motif in motif_freq:
                        motif_freq[motif] += 1
                    else:
                        motif_freq[motif] = 1
                else:
                    print('motif not valid: ', motif)

    return motif_freq
