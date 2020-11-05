"""Unit tests for infer_single_paired.py"""

from htsinfer.infer_count_motifs import count_motifs
import pytest

test_sequences_DNA = [
    "CGTACTTGGTTCGATCAACCAGGCCGTAAAGTACGACGTCATCAAAAACGTTTAAACAAA",
    "GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC",
    "GCTGCTGGTCTTCATCCAACCTATGCTCGAACCATCGGTATTTCAGTTGATCATCGACGA",
    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"
    ]

test_sequences_RNA = [
    "CGUACUUGGUUCGAUCAACCAGGCCGUAAAGUACGACGUCAUCAAAAACGUUUAAACAAA",
    "GCUCAACGUGUUGCUCCUCGACCAGCCAAAGGUUCAUUACGGCCAGUUGUUCGUGGUACC",
    "GCUGCUGGUCUUCAUCCAACCUAUGCUCGAACCAUCGGUAUUUCAGUUGAUCAUCGACGA",
    "ACACGUCGUUAUAACAUGAAAGUACGUUCUGGACGCGGUUUUUCUUUGGAUGAAAUUCGU"
    ]

test_sequences_invalid = [
    "CGTAXTTGGTTCGATRAACCAGGCCGTAPAAGTACGACMTCATQCAAAAACGTTNAAACAAA",
    "GCTGCTGGTCTTCATXCCAACCTATGCTCGAACCATCGGTATTTCAGTTGATWCATCGACGA"
    ]

test_sequences_one_wrong = [
    "CGTACTTGGTTCGATCAACCAGGCCGTAAAGTACGACGTCATCAAAAACGTTTAAACAAA",
    "GCTCAACGTGTTGCTCCTCXGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC",
    "GCTGCTGGTCTTCATCCAACCTATGCTCGAACCATCGGTATTTCAGTTGATCATCGACGA",
    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"
    ]

test_not_list = 'GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC'

test_empty_list = []

test_empty_string = ""


class TestInfer:
    """Test `count_motifs()` function."""

    def test_dna(self):
        """Test DNA sequences"""
        motifs_dict = count_motifs(
            test_sequences_DNA, min_motif_length=1, max_motif_length=6,
            nucleic_acid='dna')
        assert bool(motifs_dict) is True

    def test_rna(self):
        """Test RNA sequences"""
        motifs_dict = count_motifs(
            test_sequences_RNA, min_motif_length=1, max_motif_length=6,
            nucleic_acid='rna')
        assert bool(motifs_dict) is True

    def test_invalid_sequence_include(self):
        """Test invalid sequences and include invalids"""
        motifs_dict = count_motifs(
                test_sequences_invalid, min_motif_length=3, max_motif_length=6,
                nucleic_acid='dna', non_nucleotide_characters='include')
        assert bool(motifs_dict) is True

    def test_invalid_sequence_ignore_char(self):
        """Test invalid sequences and ignore invalid chars"""
        motifs_dict = count_motifs(
                test_sequences_invalid, min_motif_length=3, max_motif_length=6,
                nucleic_acid='dna', non_nucleotide_characters='ignore_chars')
        assert bool(motifs_dict) is True

    def test_invalid_sequence_ignore_seq(self):
        """Test invalid sequences and ignore invalid seqs"""
        motifs_dict = count_motifs(
                test_sequences_invalid, min_motif_length=3, max_motif_length=6,
                nucleic_acid='dna', non_nucleotide_characters='ignore_seqs')
        assert bool(motifs_dict) is False

    def test_partial_invalid_sequence(self):
        """Test sequence where only one is wrong and will be skipped"""
        motifs_dict = count_motifs(
            test_sequences_one_wrong, min_motif_length=3, max_motif_length=6,
            nucleic_acid='dna', non_nucleotide_characters='ignore_seqs')
        assert bool(motifs_dict) is True

    def test_not_list(self):
        """Test sequence that is not a list"""
        with pytest.raises(TypeError):
            assert count_motifs(
                test_not_list, min_motif_length=3, max_motif_length=9,
                nucleic_acid='rna')

    def test_empty_list(self):
        """Test sequence list that is empty"""
        with pytest.raises(TypeError):
            assert count_motifs(
                test_not_list, min_motif_length=3, max_motif_length=9,
                nucleic_acid='rna')

    def test_empty_string(self):
        """Test sequence that is an empty string"""
        with pytest.raises(ValueError):
            assert count_motifs(
                test_empty_list, min_motif_length=3, max_motif_length=9,
                nucleic_acid='rna')

    def test_negative_motif_length(self):
        """Test sequence with negative motif length"""
        with pytest.raises(ValueError):
            assert count_motifs(
                test_sequences_DNA, min_motif_length=-3, max_motif_length=9,
                nucleic_acid='dna')

        with pytest.raises(ValueError):
            assert count_motifs(
                test_sequences_DNA, min_motif_length=3, max_motif_length=-9,
                nucleic_acid='dna')

    def test_motif_length(self):
        """Test for min motif length > max motif length"""
        with pytest.raises(ValueError):
            assert count_motifs(
                test_sequences_DNA, min_motif_length=5, max_motif_length=2,
                nucleic_acid='dna')

    def test_min_motif_longer_seq(self):
        """Test for min motif length is longer than shortest sequence"""
        with pytest.raises(ValueError):
            assert count_motifs(
                test_sequences_DNA, min_motif_length=66, max_motif_length=70,
                nucleic_acid='dna')

    def test_incorrect_nucleic_acid(self):
        """Test if nucleic_acid is specified incorrectly"""
        with pytest.raises(ValueError):
            assert count_motifs(
                test_sequences_DNA, min_motif_length=3, max_motif_length=7,
                nucleic_acid='test')
