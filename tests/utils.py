"""Test fixtures and utilities."""

from pathlib import Path

from htsinfer.exceptions import MetadataWarning

# test parameters
TEST_FILES_DIR = Path(__file__).resolve().parent / "files"
FILE_ADAPTER = TEST_FILES_DIR / "adapters.txt"
FILE_DUMMY = Path(".")
FILE_EMPTY = TEST_FILES_DIR / "empty.fastq"
FILE_FASTA = TEST_FILES_DIR / "single.fasta"
FILE_GZIPPED = TEST_FILES_DIR / "mixed_mates_compressed.fastq.gz"
FILE_INCONSISTENT_IDS_MIXED_UNKNOWN = (
    TEST_FILES_DIR / "inconsistent_mixed_unknown.fastq"
)
FILE_INCONSISTENT_IDS_SINGLE_MATE = (
    TEST_FILES_DIR / "inconsistent_single_mate.fastq"
)
FILE_INCONSISTENT_IDS_SINGLE_OLD_NEW = (
    TEST_FILES_DIR / "inconsistent_single_old_new.fastq"
)
FILE_INVALID_SEQ_1 = TEST_FILES_DIR / "invalid_seq_1.fastq"
FILE_INVALID_SEQ_2 = TEST_FILES_DIR / "invalid_seq_2.fastq"
FILE_MATE_1 = TEST_FILES_DIR / "first_mate.fastq"
FILE_MATE_2 = TEST_FILES_DIR / "second_mate.fastq"
FILE_REAL_SAMPLE = TEST_FILES_DIR / "real_sample.fastq"
FILE_SINGLE = TEST_FILES_DIR / "single.fastq"
FILE_SRA_SAMPLE_1 = TEST_FILES_DIR / "sra_sample_1.fastq"
FILE_SRA_SAMPLE_2 = TEST_FILES_DIR / "sra_sample_2.fastq"
FILE_UNKNOWN_SEQ_ID = TEST_FILES_DIR / "unknown_seq_id.fastq"
PACKAGE_DIR = Path(__file__).resolve().parents[1] / "htsinfer"
SEQ_ID_DUMMY = ""
SEQ_ID_MATE_1 = "@SRR11971718:6:73:941:1973#0/1"
SEQ_ID_MATE_2 = "@SRR11971718:6:73:941:1973#0/2"
SEQ_ID_SINGLE = "@SRR11971718:6:73:941:1973#0"


# helper classes
class RaiseMetadataWarning:
    """Raise ``MetadataWarning``."""
    def __init__(self, *args, **kwargs):
        raise MetadataWarning


class RaiseOSError:
    """Raise ``OSError``."""
    def __init__(self, *args, **kwargs):
        raise OSError


class RaiseValueError:
    def __init__(self, *args, **kwargs):
        raise ValueError
