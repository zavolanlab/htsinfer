import numpy as np

def randomize_nucleotide_sequence(input_sequence: str,
    number_random_seq: int = 1) -> list:
    """
    Returns randomised sequence adhering to the dinucleotide probabilites.

    Args:
        sequence: input sequence.
        number_random_seq: Number of randomised sequences returned.

    Returns:
        number_random_seq x randomised sequences of equal length to the input
        sequence in a list.

    Raises:
        ValueError: input_string contains characters other than ATGC.
    """
    pass


def make_markov_matrix(input_sequence: str) -> np.matrix:
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

def create_random_sequence(input_sequence: str,
    markov_matrix: np.matrix, number_random_seq: int = 1) -> list:
    """
    Returns a list with number_random_seq randomized strings, each with length
        of the input string, and based on the probabilities from
        the markov_matrix.

    Args:
        input_sequence: Input sequence.
        markov_matrix: Probability matrix.
        number_random_seq: Number of randomised sequences returned.

    Returns:
        List of randomised sequences with number_random_seq entries.
    """
    pass

"""
Comments:
* What about non ATGC letters in the input sequence? Should the matrix be 
dynamically generated, matching other letters in the input sequence?
"""