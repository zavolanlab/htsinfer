from collections import Counter
from typing import List

import numpy as np


def randomize_nucleotide_sequence(
        input_sequences: List,
        number_random_seq: int = 1,
        min_prob: float = 0) -> List:
    """
    Returns randomised sequences adhering to the dinucleotide probabilites
    found in the input sequences.

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
        ValueError: The percentage of N's one of the input sequences is too
        high
    """

    # Error if input sequences list is empty
    if len(input_sequences) == 0:
        raise ValueError("List of input sequences is empty")

    input_seqs_upper = []
    for i, sequence in enumerate(input_sequences):
        # Error if any of the input is not a str
        if not isinstance(sequence, str):
            raise ValueError(
                f"Input sequence {i} is not of type 'str'\n"
                f"{sequence}"
            )

        # Convert input strings to upper case
        sequence = sequence.upper()
        input_seqs_upper.append(sequence)

        # Error if empty sequence is entered
        if len(sequence) == 0:
            raise ValueError(
                f"Sequence {i} is empty"
            )
        # Count the occurences of N's
        count = Counter(sequence)
        if count["N"] >= 0.05 * len(sequence):
            raise ValueError(
                f"Percentage of N in input sequence"
                f"{input_sequences.index(sequence)} is "
                f"over the limit of 5 percent. Sequence: \n"
                f"{sequence}"
                )

        # Check if there are chars other than "ATGCN"
        for i, char in enumerate(sequence):
            if char not in "ATGCN":
                raise ValueError(f"Input sequence "
                                 f"{input_sequences.index(sequence)}"
                                 f"contains characters other "
                                 f"than 'ATGCN'.\n"
                                 f"First occurence:\n"
                                 f"{sequence} \n"
                                 f'{"^":>{i+1}}'
                                 )

    input_sequences = input_seqs_upper

    markov_matrix = make_markov_matrix(input_sequences)

    random_sequences = create_random_sequence(
        input_seqs=input_sequences,
        markov_matrix=markov_matrix,
        number_random_seq=number_random_seq)

    return random_sequences


def make_markov_matrix(sequences: list):
    """
    Returns markov matrix based on the input sequence,

          A      C      G      T
        A P(AA)  P(AC)  P(AG)  P(AT)
        C P(TA)  P(CC)  P(TG)  P(CT)
        G P(GA)  P(GC)  P(GG)  P(GT)
        T P(CA)  P(TC)  P(CG)  P(TT)

        where P(AA) is the probability of A following an A.

    Args:
        input_sequence: Input sequence.

    Returns:
        Markov matrix as 4x4 np.matrix with propabilities of nucleotide pairs.

    Raises:
        ValueError: input_string contains characters other than ATGC.
    """
    nucl = ["A", "C", "G", "T"]

    # Possible nucleotide pairs
    comb = [["AA", "AC", "AG", "AT"],
            ["CA", "CC", "CG", "CT"],
            ["GA", "GC", "GG", "GT"],
            ["TA", "TC", "TG", "TT"]]
    result = np.ones((4, 4))

    # count occurence of nucleotide pairs:
    for seq in sequences:
        for i in range(len(nucl)):  # i represents the row
            for j in range(len(nucl)):  # j represents the column
                a = CountOccurrences(seq,  comb[i][j])
                result[i, j] += a
                a = 1

    # convert occurence to probabilities
    for i in range(len(nucl)):  # i represents the row
        tot = np.sum(result[i])
        for j in range(len(nucl)):  # j represents the column
            if result[i, j] != 0:
                prob = result[i, j]/tot
                result[i, j] = prob

    return result


def CountOccurrences(string, substring):

    # Initialize count and start to 0
    count = 0
    start = 0

    # Search through the string till
    # we reach the end of it
    while start < len(string):

        # Check if a substring is present from
        # 'start' position till the end
        pos = string.find(substring, start)

        if pos != -1:
            # If a substring is present, move 'start' to
            # the next position from start of the substring
            start = pos + 1

            # Increment the count
            count += 1
        else:
            # If no further substring is present
            break
    # return the value of count
    return count


def create_random_sequence(
        input_seqs: List,
        markov_matrix: np.array,
        number_random_seq: int = 1) -> List:
    """
    Returns a list with number_random_seq randomised strings for each input
        string in the input_seqs list, each with length of the corresponding
        input string and based on the probablilities from the markov_matrix.
        The key is used to link the matrix values to the bases: e.g. the matrix
          A      C      G      T
        A P(AA)  P(AC)  P(AG)  P(AT)
        C P(TA)  P(CC)  P(TG)  P(CT)
        G P(GA)  P(GC)  P(GG)  P(GT)
        T P(CA)  P(TC)  P(CG)  P(TT)
        has the corresponding key ("A","C","G","T").

    Args:
        input_seqs: List of input sequences.
        markov_matrix: Markov matrix of the probabilities of nucleotide pairs.
        number_random_seq: Number of random sequences to return per input seq.

    Returns:
        List of lists with random sequences seuqences, each list is based on
        one input string. The length of each list corresponds to
        number_random_seq.
    """
    # Key to read the markov matrix and chose random nucleotides
    key = ("A", "C", "G", "T")
    random_list_all = []

    # Iterate through all input_seqs
    for input_string in input_seqs:
        random_list = []

        # For each input_seq create number_random_seq random sequences
        for _ in range(number_random_seq):
            random_string = np.random.choice(a=key)
            for i in range(len(input_string)-1):
                # Add the remaining letter based on the markov matrix.
                # key is used as pool of possible choices and to find the
                # row of the matrix corresponding to the correct base.
                random_string += np.random.choice(
                    a=key,
                    p=markov_matrix[key.index(random_string[i])])
            random_list.append(random_string)
        random_list_all.append(random_list)

    return(random_list_all)


def Markov_check(Generated_array, Seq_array, treshold):
    """
    Goal: Compare the markov matrix of the input and output sequences.

    Currently this is not implemented in the main function.
    """
    # Input: arrays with nucleotide frequences in a 4x4 array
    # Output: String giving the nucleotide pairs with a frequency
    # deviation higher than the specified treshold
    # Array shape
    comb = [["AA", "AC", "AG", "AT"],
            ["CA", "CC", "CG", "CT"],
            ["GA", "GC", "GG", "GT"],
            ["TA", "TC", "TG", "TT"]]
    # Used nucleotides
    nucl = ["A", "C", "G", "T"]
    # Output error list
    Errorlist = []
    # Error treshold
    treshold = treshold
    # treshold_per = round((treshold - 1) * 100,1)

    # Division of original markov array (seq array) used to generate
    # the sequence (generated array) to calculate freq difference
    Dev_array = np.divide(Generated_array, Seq_array)

    # checking for more than 5% variation in dinucleotide occurences
    for i in range(len(nucl)):  # i represents the row 
        for j in range(len(nucl)):  # j represents the column
            if Dev_array[i][j] >= treshold-1:
                Errorlist.append(comb[i][j])
            else:
                pass
                # d = 0

    if not Errorlist:
        return True
    else:
        return False
    """
    if not Errorlist :
        Markov_check = f'There is not more than a {treshold_per}% deviation in the frequences of the computed nucleotide string with comparison to the original dataset'
    else:
        Markov_check = f'There is more than a {treshold_per}% deviation in the frequences of the following nucleotides: {Errorlist}' 
    
    return Markov_check
    """

"""
Limiting case: What to do when certain di-nucleotide combinations never occur.
"""
