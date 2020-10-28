import numpy as np
from collections import Counter

def randomize_nucleotide_sequence(input_sequences: list,
    number_random_seq: int = 1, min_prob: float = 0) -> list:
    """
    Returns randomised sequence adhering to the dinucleotide probabilites.

    Args:
        input_sequences: List of input sequences as strings.
        number_random_seq: Number of randomised sequences returned per input
            sequence.
        min_prob: Minimal probability of dinucleotide to appear. If set to -1
            probability of 0 is possible. If set to 0, minimum probability will
            correspond to one count in all input sequences.

    Returns:
        number_random_seq x randomised sequences of equal length to the input
        sequence in a list.

    Raises:
        ValueError: input_string contains too many characters other than ATGCN.
        ValueError: The percentage of N's one of the input sequences is too high
    """
    
    for sequence in input_sequences:
        # Count the occurences of N's
        count = Counter(sequence)
        if count["N"] >= 0.05 * len(sequence):
            raise ValueError(
                f"Percentage of N in input sequence  {input_sequences.index(sequence)} is "
                f"over the limit of 5 percent. Sequence: \n"
                f"{sequence}"
        )

        # Check if there are chars other than "ATGCN"
        for i,char in enumerate(sequence):
            if char not in "ATGCN":
                correct = 'correct'
                raise ValueError(f"Input sequence "
                f"{input_sequences.index(sequence)} contains characters other "
                f"than 'ATGCN'.\n"
                f"First occurance:\n"
                f"{sequence} \n"
                f'{"^":>{i+1}}'
            )


def make_markov_matrix(input_sequence: list) -> np.array:
    """
    Returns markov matrix based on the input sequence,

          A      T      G      C
        A P(AA)  P(AT)  P(AG)  P(AC)
        T P(TA)  P(TT)  P(TG)  P(TC)
        G P(GA)  P(GT)  P(GG)  P(GC)
        C P(CA)  P(CT)  P(CG)  P(CC)

        where P(AA) is the probability of A following an A.

    Args:
        input_sequence: Input sequence.

    Returns:
        Markov matrix as 4x4 np.matrix with propabilities of nucleotide pairs.

    Raises:
        ValueError: input_string contains characters other than ATGC.
    """
    pass


def create_random_sequence(input_seqs: list, markov_matrix: np.array,
    number_random_seq: int = 1) -> list:
    """
    Returns a list with number_random_seq randomised strings for each input 
        string in the input_seqs list, each with length of the corresponding
        input string and based on the probablilities from the markov_matrix.
        The key is used to link the matrix values to the bases: e.g. the matrix
                A      T      G      C
            A P(AA)  P(AT)  P(AG)  P(AC)
            T P(TA)  P(TT)  P(TG)  P(TC)
            G P(GA)  P(GT)  P(GG)  P(GC)
            C P(CA)  P(CT)  P(CG)  P(CC)
        has the corresponding key ("A","T","G","C").

    Args:
        input_seqs: List of input sequences.
        markov_matrix: Markov matrix of the probabilities of nucleotide pairs.
        number_random_seq: Number of random sequences return per input sequence.

    Returns:
        List of lists with random sequences seuqences, each list is based on
        one input string. The length of each list corresponds to
        number_random_seq.
    """
    # Key to read the markov matrix and chose random nucleotides
    key = ("A", "T", "G", "C")
    random_list_all = []

    # Iterate through all input_seqs
    for input_string in input_seqs:
        random_list = []

        # For each input_seq create number_random_seq random sequences
        for _ in range(number_random_seq):
            random_string = np.random.choice(key) #first char is random from key
            for i in range(len(input_string)-1):
                # Add the remaining letter based on the markov matrix.
                # key is used as pool of possible choices and to find the 
                # row of the matrix corresponding to the correct base.
                random_string += np.random.choice(a = key, 
                    p = markov_matrix[key.index(random_string[i])])
            random_list.append(random_string)
        random_list_all.append(random_list)
    
    return(random_list_all)




"""
Comments:
* Add some sort control to check, whether the random sequences produce the 
same markov matrix as the input sequences.
* Maybe we need to deal with small letters.

Mischa to-do:
* Implement ValueError: Too many N's.

Thomas to-do:
* Update the markov-generation, make sure it works with the N's and stuff.
* Implement function to compare the markov matrix of the in and output.

Both:
* Start with implementing tests.
"""