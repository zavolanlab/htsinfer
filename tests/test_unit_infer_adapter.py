"""Unit tests for infer_adapter.py"""

import os
from htsinfer.infer_adapter import infer

path = os.path.dirname(__file__)
path = os.path.join(path, "test_files")


def test_single_file():
    """Test single file"""
    file1 = os.path.join(path, "first_mate.fastq")
    result_1, result_2 = infer(file1)
    assert result_1 == "NA" and \
        result_2 == "not_available"


def test_compressed_file():
    """Test single compressed file."""
    file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
    result_1, result_2 = infer(file1)
    assert result_1 == "NA" and \
        result_2 == "not_available"


def test_empty_file():
    """Test empty file."""
    file1 = os.path.join(path, "SRR11971713_1.fastq")
    file2 = os.path.join(path, "empty.fastq")
    result_1, result_2 = infer(file_1=file1, file_2=file2)
    assert result_1 == "NA" and \
        result_2 == "NA"


def test_adapter():
    """Test adapter name."""
    file1 = os.path.join(path, "SRR13492831.fastq")
    result_1, result_2 = infer(file_1=file1)
    assert result_1 == "AAAAAAAAAAAAAAA" and \
        result_2 == "not_available"


def test_max_records():
    """Set max_records argument"""
    file1 = os.path.join(path, "SRR13492831.fastq")
    result_1, result_2 = infer(file_1=file1, max_records=2)
    assert result_1 == "NA" and \
        result_2 == "not_available"


def test_confidence():
    """Test confidence score"""
    file1 = os.path.join(path, "SRR13492831.fastq")
    result_1, result_2 = infer(file_1=file1, min_match=20, factor=3)
    assert result_1 == "NA" and \
        result_2 == "not_available"


def test_invalid_file():
    """Test invalid file"""
    file1 = os.path.join(path, "SRR11972514_1.file_1")
    result_1, result_2 = infer(file_1=file1)
    assert result_1 == "invalid_file" and \
        result_2 == "not_available"


def test_min_match():
    """Test min match percentage"""
    file1 = os.path.join(path, "SRR13496438.fastq.gz")
    result_1, result_2 = infer(file_1=file1, min_match=70)
    assert result_1 == "NA" and \
        result_2 == "not_available"


def test_params():
    """Test all params"""
    file1 = os.path.join(path, "SRR13496438.fastq.gz")
    file2 = os.path.join(path, "SRR13492831.fastq")
    result_1, result_2 = infer(
        file_1=file1,
        file_2=file2,
        min_match=30,
        factor=1.2
        )
    assert result_1 == "GATCGGAAGAGCACA" and \
        result_2 == "AAAAAAAAAAAAAAA"
