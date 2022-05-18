"""Test fixtures and utilities."""

from pathlib import Path

from htsinfer.exceptions import (FileProblem, MetadataWarning)
from htsinfer.models import Source

# test files
PACKAGE_DIR = Path(__file__).resolve().parents[1] / "htsinfer"
TEST_FILES_DIR = Path(__file__).resolve().parent / "files"
FILE_ADAPTER = TEST_FILES_DIR / "adapter_fragments.txt"
FILE_ADAPTER_INVALID_CHARS = (
    TEST_FILES_DIR / "adapter_fragments_invalid_chars.txt"
)
FILE_ADAPTER_SEQ_TOO_LONG = (
    TEST_FILES_DIR / "adapter_fragments_seq_too_long.txt"
)
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
FILE_ORIENTATION_ISF_1 = TEST_FILES_DIR / "orientation_isf_mate_1.fastq"
FILE_ORIENTATION_ISF_2 = TEST_FILES_DIR / "orientation_isf_mate_2.fastq"
FILE_ORIENTATION_ISR_1 = TEST_FILES_DIR / "orientation_isr_mate_1.fastq"
FILE_ORIENTATION_ISR_2 = TEST_FILES_DIR / "orientation_isr_mate_2.fastq"
FILE_ORIENTATION_IU_1 = TEST_FILES_DIR / "orientation_iu_mate_1.fastq"
FILE_ORIENTATION_IU_2 = TEST_FILES_DIR / "orientation_iu_mate_2.fastq"
FILE_ORIENTATION_SF = TEST_FILES_DIR / "orientation_sf.fastq"
FILE_ORIENTATION_SR = TEST_FILES_DIR / "orientation_sr.fastq"
FILE_ORIENTATION_U = TEST_FILES_DIR / "orientation_u.fastq"
FILE_2000_RECORDS = TEST_FILES_DIR / "2000_records.fastq"
FILE_SINGLE = TEST_FILES_DIR / "single.fastq"
FILE_SRA_SAMPLE_1 = TEST_FILES_DIR / "sra_sample_1.fastq"
FILE_SRA_SAMPLE_2 = TEST_FILES_DIR / "sra_sample_2.fastq"
FILE_TRANSCRIPTS = TEST_FILES_DIR / "transcripts.fa"
FILE_TRANSCRIPTS_GZ = TEST_FILES_DIR / "transcripts.fa.gz"
FILE_UNKNOWN_SEQ_ID = TEST_FILES_DIR / "unknown_seq_id.fastq"
FILE_UNMAPPED_PAIRED_1 = TEST_FILES_DIR / "unmapped_paired_mate_1.fastq"
FILE_UNMAPPED_PAIRED_2 = TEST_FILES_DIR / "unmapped_paired_mate_2.fastq"
FILE_UNMAPPED_SINGLE = TEST_FILES_DIR / "unmapped_single.fastq"

# test parameters
DICT_DF = {
    'key_3': 1,
    'key_1': 2,
    'key_2': 3,
}
SOURCE_HUMAN = Source(short_name="hsapiens", taxon_id=9606)
SEQ_ID_DUMMY = ""
SEQ_ID_MATE_1 = "@SRR11971718:6:73:941:1973#0/1"
SEQ_ID_MATE_2 = "@SRR11971718:6:73:941:1973#0/2"
SEQ_ID_SINGLE = "@SRR11971718:6:73:941:1973#0"


# helper classes
class RaiseFileProblem:
    def __init__(self, *args, **kwargs):
        raise FileProblem


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
