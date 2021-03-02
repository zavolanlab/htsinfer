"""Unit tests for infer_organism.py"""

import os
from htsinfer.infer_organism import infer

dir_path = os.path.abspath(os.path.join(__file__, "../.."))
path = os.path.dirname(__file__)
transcript_fasta_path = os.path.join(
    dir_path, "htsinfer", "transcripts.fasta.zip"
    )
path = os.path.join(path, "test_files")


def test_single_file():
    """Test single file"""
    file1 = os.path.join(path, "first_mate.fastq")
    result = infer(transcript_fasta=transcript_fasta_path, file_1=file1)
    assert result == "invalid_file"


def test_compressed_file():
    """Test single compressed file."""
    file1 = os.path.join(path, "SRR13496438.fastq.gz")
    result = infer(transcript_fasta=transcript_fasta_path, file_1=file1)
    assert result == "NA"


def test_empty_file():
    """Test empty file."""
    file1 = os.path.join(path, "empty.fastq")
    result = infer(transcript_fasta=transcript_fasta_path, file_1=file1)
    assert result == "invalid_file"


def test_organism():
    """Test organism name."""
    file1 = os.path.join(path, "SRR13496438.fastq.gz")
    result = infer(
        transcript_fasta=transcript_fasta_path, file_1=file1, factor=1
        )
    assert result == "oaries"


def test_min_match():
    """Test min match percentage"""
    file1 = os.path.join(path, "SRR13496438.fastq.gz")
    result = infer(
        transcript_fasta=transcript_fasta_path, file_1=file1, min_match=70
        )
    assert result == "NA"


def test_params():
    """Test all params"""
    file1 = os.path.join(path, "SRR13496438.fastq.gz")
    file2 = os.path.join(path, "SRR13492831.fastq")
    result = infer(
        transcript_fasta=transcript_fasta_path,
        file_1=file1,
        file_2=file2,
        min_match=10,
        factor=1
        )
    assert result == "oaries"
