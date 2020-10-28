"""Unit tests for extract_random_seq_subset.py"""

import os
import stat
import filecmp
import shutil
from htsinfer.extract_random_seq_subset import extract

path = os.path.dirname(__file__)
path = os.path.join(path, "test_files")
outpath = os.path.join(path, "outputs")


class TestExtract:
    """Test `extract()` function."""

    # set up output directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    def test_no_input_output_file(self):
        """Test no input and no output file provided"""
        assert extract([], []) == "no_input_output_provided"

    def test_no_input_file(self):
        """Test no input file provided"""
        file1 = os.path.join(outpath, "SRR11971558_1_subset.fastq")
        assert extract([], [file1]) == "no_input_provided"

    def test_no_output_file(self):
        """Test no output file provided"""
        file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
        assert extract([file1], []) == "no_output_provided"

    def test_output_file_error(self):
        """Test unwritable output file"""
        file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
        outfile1 = os.path.join(outpath, "no_seq.fastq")
        os.chmod(outfile1, stat.S_IREAD)
        outcome = extract([file1], [outfile1])
        os.remove(outfile1)
        assert outcome == "file_error"

    def test_mismatched_in_out_file(self):
        """Test mismatched number of  in/output files"""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        file2 = os.path.join(path, "SRR11971713_1.fastq")
        outfile1 = os.path.join(outpath, "SRR11971558_2_subset.fastq")
        assert extract([file1, file2], [outfile1]) == "mismatched_in_out_files"

    def test_same_in_out_files(self):
        """Input and output files are the same"""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        file2 = os.path.join(outpath, "SRR11971558_2.fastq")
        shutil.copy2(file1, file2)
        assert extract([file2], [file2]) == "same_in_out_files"

    def test_no_seq_found(self):
        """Test no input sequences found"""
        file1 = os.path.join(path, "no_seq.fastq")
        outfile1 = os.path.join(outpath, "no_seq.fastq")
        assert extract([file1], [outfile1]) == "no_input_sequences"

    def test_empty_file(self):
        """Test empty file"""
        file1 = os.path.join(path, "empty.fastq")
        outfile1 = os.path.join(outpath, "empty.fastq")
        outcome = extract([file1], [outfile1])
        os.remove(outfile1)
        assert outcome == "no_input_sequences"

    def test_file_not_found(self):
        """Test file not found"""
        file1 = os.path.join(path, "missing.fastq")
        outfile1 = os.path.join(outpath, "missing.extracted.fastq")
        assert extract([file1], [outfile1]) == "file_error"

    def test_compressed_file(self):
        """Test single compressed file."""
        file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
        outfile1 = os.path.join(outpath, "subset_from_compressed.fastq")
        assert extract([file1], [outfile1]) == "extraction_done"

    def test_copying_file(self):
        """Test copy input file."""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath, "subset_all.fastq")
        extract([file1], [outfile1], max_records=0)
        assert filecmp.cmp(file1, outfile1, shallow=False)

    def test_extract_proportion(self):
        """Test extraction of a proportion of reads."""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath, "subset_proportion.fastq")

        prop = 0.4
        outcome = extract([file1], [outfile1], proportion=prop)

        (rowin, rowout) = (0, 0)

        with open(file1) as f:
            rowin = len(f.readlines())
        with open(outfile1) as f:
            rowout = len(f.readlines())

        assert outcome == "extraction_done" and int(rowin * prop) == rowout

    def test_extract_too_large_proportion(self):
        """Test extraction of a proportion of reads."""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath, "subset_proportion.fastq")

        prop = 2.0
        outcome = extract([file1], [outfile1], proportion=prop)

        (rowin, rowout) = (0, 0)

        with open(file1) as f:
            rowin = len(f.readlines())
        with open(outfile1) as f:
            rowout = len(f.readlines())

        assert outcome == "extraction_done" and rowin == rowout

    def test_extract_number(self):
        """Test extraction of a certain number of records."""
        lines_per_rec = 4  # fastq format

        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath,
                                "SRR11971558_2_subset_OK_number.fastq")

        num = 1

        outcome = extract([file1], [outfile1], max_records=num)

        rowout = 0

        with open(outfile1) as f:
            rowout = len(f.readlines())

        assert (outcome == "extraction_done"
                and int(rowout/lines_per_rec) == num)

    def test_extract_too_large_number(self):
        """Test extraction when requested records
             is larger than available records."""

        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath,
                                "SRR11971558_2_subset_too_large_number.fastq")

        num = 100

        outcome = extract([file1], [outfile1], max_records=num)

        (rowin, rowout) = (0, 0)

        with open(file1) as f:
            rowin = len(f.readlines())
        with open(outfile1) as f:
            rowout = len(f.readlines())

        assert outcome == "extraction_done" and rowin == rowout

    def test_invalid_max_records(self):
        """Test invalid numerical parameters."""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath, "SRR11971558_2_max_records-1.fastq")
        assert (extract([file1], [outfile1], max_records=-1)
                == "invalid_number_requested")

    def test_invalid_proportion(self):
        """Test invalid numerical parameters."""
        file1 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath, "SRR11971558_2_proportion-1.fastq")
        assert (extract([file1], [outfile1], proportion=-1.0)
                == "invalid_number_requested")

    def test_multiple_files(self):
        """Test multiple files."""
        file1 = os.path.join(path, "SRR11971558_1.fastq.gz")
        file2 = os.path.join(path, "SRR11971558_2.fastq")
        outfile1 = os.path.join(outpath, "SRR11971558_1_M1.fastq")
        outfile2 = os.path.join(outpath, "SRR11971558_2_M1.fastq")
        assert (extract([file1, file2], [outfile1, outfile2])
                == "extraction_done")
