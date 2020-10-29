"""Test for rand_nuc_seq.py"""

import pytest
from htsinfer.rand_nuc_seq import randomize_nucleotide_sequence, make_markov_matrix


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
        test_output = randomize_nucleotide_sequence(test_input, number_random_seq=10)
        combined = []
        for seqs in test_output:
            combined.extend(seqs)
        assert 30 == len(combined)

'''
    def test_markov_matrix_probabilities(self):
        """Checks if the markov matrix probabilities are correct"""
        test_input = ["AT","TT","GT","CT"]
        markov_matrix = make_markov_matrix(test_input)

'''