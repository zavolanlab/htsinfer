import numpy as np

def randomize_nucleotide_sequence(input_sequence: str) -> str:
    """
    Returns randomised sequence adhering to the dinucleotide probabilites.

    Args:
        sequence: input sequence.

    Returns:
        Randomised sequence of equal length to the input sequence.
    """


def make_markov_matrix(input_sequence: str) -> np.matrix:
    """
    Returns markov matrix based on the input sequence.
      A      T   G C
    A P(AA)  AT  
    T ... 
    G
    C

    Args:
        input_sequence: Input sequence.

    Returns:
        Markov matrix as np.matrix with propabilities of nucleotide pairs.

    
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

"""
Comments:
* What about non ATGC letters in the input sequence?

"""