"""Unit tests for infer_single_paired.py"""

import os
from htsinfer.infer_single_paired import infer

path = os.path.dirname(__file__)
path = os.path.join(path, "sample_files")


def test_single_file():
    file_1 = os.path.join(path, "first_mate.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file_1)
    assert result_lib_1 == "first_mate" and \
        result_lib_2 == "not_available" and \
        mate_relationship == "not_available"


def test_compressed_file():
    file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
    result_lib_1, result_lib_2, mate_relationship = infer(file1)
    assert result_lib_1 == "mixed_mates" and \
        result_lib_2 == "not_available" and \
        mate_relationship == "not_available"


def test_second_mate():
    file1 = os.path.join(path, "SRR11971713_2.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file1)
    assert result_lib_1 == "second_mate" and \
        result_lib_2 == "not_available" and \
        mate_relationship == "not_available"


def test_split_mates():
    file1 = os.path.join(path, "first_mate.fastq")
    file2 = os.path.join(path, "second_mate.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file1, file2)
    assert result_lib_1 == "first_mate" and \
        result_lib_2 == "second_mate" and \
        mate_relationship == "split_mates"


def test_invalid_identifiers():
    file1 = os.path.join(path, "SRR11971713_1.fastq")
    file2 = os.path.join(path, "SRR11971558_2.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file1, file2)
    assert result_lib_1 == "no_mate_info" and \
        result_lib_2 == "no_mate_info" and \
        mate_relationship == "not_mates"


def test_empty_file():
    file1 = os.path.join(path, "empty.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file1)
    assert result_lib_1 == "invalid_file" and \
        result_lib_2 == "not_available" and \
        mate_relationship == "not_available"


def test_not_mates():
    file1 = os.path.join(path, "SRR11972514_1.fastq")
    file2 = os.path.join(path, "SRR11971713_2.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file1, file2)
    assert result_lib_1 == "first_mate" and \
        result_lib_2 == "second_mate" and \
        mate_relationship == "not_mates"


def test_invalid_file():
    file1 = os.path.join(path, "SRR11972514_1.fast")
    result_lib_1, result_lib_2, mate_relationship = infer(file1)
    assert result_lib_1 == "invalid_file" and \
        result_lib_2 == "not_available" and \
        mate_relationship == "not_available"


def test_max_records():
    file1 = os.path.join(path, "first_mate.fastq")
    file2 = os.path.join(path, "second_mate.fastq")
    result_lib_1, result_lib_2, mate_relationship = infer(file1, file2, 2)
    assert result_lib_1 == "first_mate" and \
        result_lib_2 == "second_mate" and \
        mate_relationship == "split_mates"
