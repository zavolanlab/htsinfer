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


def create_random_sequence(input_seqs: list, markov_matrix: np.array,
    key: str, ignore_special_char: bool = 1,
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
        has the corresponding key "ATGC".

    Args:
        input_seqs: List of input sequences.
        markov_matrix: Markov matrix of the probabilities of nucleotide pairs.
        key: Possible nucleotides and key to matrix.
        ignore_special_char: If it is set to true, in case of a non ATGC
            character in a input string, the following base will be randomly
            chosen from ATGC.
        number_random_seq: Number of random sequences return per input sequence.

    Returns:
        List of lists with random sequences seuqences, each list is based on
        one input string. The length of each list corresponds to
        number_random_seq.
    """

    



"""
Comments:
* What about non ATGC letters in the input sequence? Should the matrix be 
dynamically generated, matching other letters in the input sequence?
"""