"""Unit tests for infer_extend_motifs.py"""

import pytest
from htsinfer.infer_extend_motifs import extend_motifs

TEST_SEQUENCES_DNA = [
    "CGTACTTGGTTCGATCAACCAGGCCGTAAAGTACGACGTCGTCAAAAACGTTTAAACAAA",
    "GCTCAACGTGTTGCTCCTAGACCAGCCAAAGGTTCAATACGGCCAGTTGTTCGTGGTACC",
    "GCTGCTGGTCTTCAACCAACCTATGCTCGAACCATCGGTATTCCAATTCATCAACGACGA",
    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"
    ]

TEST_SEQUENCES_RNA = [
    "CGUACUUGGUUCGAUCAACCAGGCCGUAAAGUACGACGUCAUCAAAAACGUUUAAACAAA",
    "GCUCGACGUGUUGCUGCUCGACCAGCCAAAGUUUCAUUACGGCCAGUUGUUCGUGGUACC",
    "GCUGCUGGUCUUCAUCCAACCUAUGCUCGAACCAUCGGUAUUUCAGUUGAUCAUCGACGA",
    "ACACGUCGUUAUAACAUGAAAGUACGUUCUGGACGCGGUUUUUCUUUGGAUGAAAUUCGU"
    ]

TEST_MOTIFS_DNA = ["GCC", "GTT", "AAG", "ATT", "TCA", "TTC", "CTC", "CCA"]

TEST_MOTIFS_RNA = ["GCC", "GUU", "AAG", "AUU", "UCA", "UUC", "CUC", "CCA"]

TEST_LONG_MOTIF_DNA = ["GCCCCCCCCCCCCGGGCCCCC", "GTTATA", "AAG",
                       "ATT", "TCA", "TTC", "CTC", "CCA"]

TEST_SHORT_SEQUENCES_DNA = [
    "CGTACTTGGT",
    "GCTCAACGTGTTGCTCCTAGACCAGCCAAAG",
    "GCTGCTGGTCTTCAACCAACCTATGCTCGAACCATCGGTATTCC",
    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"
    ]

TEST_NOT_LIST = 'GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC'

TEST_EMPTY_LIST: list = []

TEST_EMPTY_STRING_MOTIF = ""

TEST_EMPTY_STRING_SEQUENCE = ""


class TestInfer:
    """Test `extend_motifs()` function."""

    def test_dna(self):
        """Test DNA sequences"""
        extended_motifs = extend_motifs(TEST_SEQUENCES_DNA,
                                        TEST_MOTIFS_DNA, 'dna')
        assert bool(extended_motifs) is True

    def test_rna(self):
        """Test RNA sequences"""
        extended_motifs = extend_motifs(TEST_SEQUENCES_RNA,
                                        TEST_MOTIFS_RNA, 'rna')
        assert bool(extended_motifs) is True

    def test_not_a_list(self):
        """Test sequence that is not a list"""
        with pytest.raises(TypeError):
            assert extend_motifs(TEST_NOT_LIST,
                                 TEST_MOTIFS_DNA, nucleic_acid='dna')

    def test_empty_list(self):
        """Test sequence list that is empty"""
        with pytest.raises(ValueError):
            assert extend_motifs(TEST_EMPTY_LIST,
                                 TEST_MOTIFS_RNA, nucleic_acid='rna')

    def test_empty_string_motifs(self):
        """Test sequence (from input_sequences) that is an empty string"""
        with pytest.raises(ValueError):
            assert extend_motifs(TEST_SEQUENCES_RNA,
                                 TEST_EMPTY_STRING_MOTIF, nucleic_acid='rna')

    def test_empty_string_sequences(self):
        """Test sequence (from input_motifs) that is an empty string"""
        with pytest.raises(ValueError):
            assert extend_motifs(TEST_EMPTY_STRING_SEQUENCE,
                                 TEST_MOTIFS_RNA, nucleic_acid='rna')

    def test_min_motif_longer_seq(self):
        """Test for min motif length is longer than shortest sequence"""
        with pytest.raises(ValueError):
            assert extend_motifs(TEST_SHORT_SEQUENCES_DNA,
                                 TEST_LONG_MOTIF_DNA, nucleic_acid='dna')

    def test_incorrect_nucleic_acid(self):
        """Test if nucleic_acid is specified incorrectly"""
        with pytest.raises(ValueError):
            assert extend_motifs(TEST_SEQUENCES_DNA,
                                 TEST_MOTIFS_DNA, nucleic_acid='test')
