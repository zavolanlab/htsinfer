"""Unit tests for infer_single_paired.py"""

import os
import pytest
import numpy as np
from htsinfer.infer_read_layout import infer  # noqa: F401
from htsinfer.infer_read_layout import (randomize_nucleotide_sequence,
                                        make_markov_matrix)

path = os.path.dirname(__file__)
path = os.path.join(path, "test_files")


class TestRandNuqSeq:
    """Tests 'randomize_nucleotide_sequence()' function"""

    def test_n_error(self):
        """Checks for ValueError if percentage of N's is too high"""
        test_input = ["ATGCATGACT", "ACGTACGN"]
        with pytest.raises(ValueError):
            randomize_nucleotide_sequence(test_input)

    def test_unaccepted_char(self):
        """Checks if chars other than ATGCN are in the input sequence"""
        test_input = ["ATGCATGATCGAT", "ACGTAGGGAT", "ATGCATGu"]
        with pytest.raises(ValueError):
            randomize_nucleotide_sequence(test_input)

    def test_output_list_length(self):
        """Checks if the returned list has the excpected length"""
        test_input = ["ATGT", "ATCGT", "A"]
        test_output = randomize_nucleotide_sequence(test_input,
                                                    number_random_seq=10)
        combined = []
        for seqs in test_output:
            combined.extend(seqs)
        assert len(combined) == 30

    def test_number_random_seq_error(self):
        """Checks ValueError if the number of random seqs is > 1"""
        test_input = ["AGATG", "AT"]
        with pytest.raises(ValueError):
            randomize_nucleotide_sequence(test_input, number_random_seq=-4)

    def test_pseudo_count_error(self):
        """Checks ValueError for negative pseudo counts"""
        test_input = ["AGAT"]
        with pytest.raises(ValueError):
            randomize_nucleotide_sequence(test_input, pseudo_count=-3)

    def test_markov_matrix_probabilities(self):
        """Checks if the markov matrix probabilities are correct"""
        test_input = ["AT", "TT", "GT", "CT"]
        markov_matrix = make_markov_matrix(test_input, pseudo_count=1)
        test_matrix = np.array([[0.2, 0.2, 0.2, 0.4],
                                [0.2, 0.2, 0.2, 0.4],
                                [0.2, 0.2, 0.2, 0.4],
                                [0.2, 0.2, 0.2, 0.4]])
        assert (markov_matrix == test_matrix).all()
