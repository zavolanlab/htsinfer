"""Unit Tests for infer_single_paired.py"""

import os
import pytest
from src.infer_single_paired import EndParser


path = os.path.dirname(__file__)
path = os.path.join(path, "sample_files")


@pytest.mark.skip()
def test_invalid_identifier():
    """ Invalid Identifier """
    file1 = os.path.join(path, "SRR11972507_1.fastq")
    p1, p2 = obj.fastq(file1)
    assert p1 == 0 and p2 == -1


@pytest.mark.skip()
def test_empty():
    """ Empty File """
    file1 = os.path.join(path, "empty.fastq")
    p1, p2 = obj.fastq(file1)
    assert p1 == 0 and p2 == -1


@pytest.mark.skip()
def test_mixed_identifier_1():
    """ Mixed - Second Mate : Identifier type 1"""
    file1 = os.path.join(path, "SRR11971718_1.fastq")
    file2 = os.path.join(path, "SRR11971713_2.fastq")
    p1, p2 = obj.fastq(file1, file2)
    assert p1 == 3 and p2 == 2


@pytest.mark.skip()
def test_single_identifier_2():
    """ Single End : Identifier type 2"""
    file1 = os.path.join(path, "SRR11972514_1.fastq")
    p1, p2 = obj.fastq(file1)
    assert p1 == 1 and p2 == -1


@pytest.mark.skip()
def test_compressed():
    """ Compressed fastq """
    file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
    p1, p2 = obj.fastq(file1)
    assert p1 == 3 and p2 == -1


@pytest.mark.skip(reason="SRA toolkit dependency")
def test_sra_files_single_double():
    """ First Mate - Second Mate """
    file1 = os.path.join(path, "sra_1.fastq")
    file2 = os.path.join(path, "sra_2.fastq")
    p1, p2 = obj.fastq(file1, file2)
    assert p1 == 1 and p2 == 2


@pytest.mark.skip(reason="SRA toolkit dependency")
def test_sra_files_single_mixed():
    """ First Mate - Mixed """
    file1 = os.path.join(path, "sra_1.fastq")
    file2 = os.path.join(path, "sra_combined.fastq")
    p1, p2 = obj.fastq(file1, file2)
    assert p1 == 1 and p2 == 3
