"""Data models."""

from enum import (
    Enum,
    IntEnum,
)
import logging
import re
from typing import Optional
from pathlib import Path
import tempfile

# pylint: disable=no-name-in-module,invalid-name
from pydantic import BaseModel


class CleanupRegimes(Enum):
    """Enumerator of cleanup regimes."""
    DEFAULT = "default"
    KEEP_ALL = "keep_all"
    KEEP_NONE = "keep_none"
    KEEP_RESULTS = "keep_results"


class LogLevels(Enum):
    """Log level enumerator."""
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARN = logging.WARNING
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL


class RunStates(IntEnum):
    """Enumerator of run states and exit codes."""
    OKAY = 0
    WARNING = 1
    ERROR = 2


SeqIdFormats = Enum(  # type: ignore
    # Source information:
    # https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    # https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
    value='SeqIdFormats',
    names=[
        # Illumina Casava >=1.8
        # example EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
        (
            'Casava >=1.8',
            re.compile(
                r'(?P<prefix>\w+:\d+:\w+:\d+:\d+:\d+:\d+(:[ACGTN]\+[ACGTN])?)'
                r'('
                r'(?P<linker> )'
                r'(?P<mate>[12])'
                r'(?P<suffix>:[YN]:\d*[02468]:([ACGTN]|\d)+)'
                r')?'
            )
        ),
        # Illumina Casava <1.8
        # example: HWUSI-EAS100R:6:73:941:1973#0/1
        (
            'Casava <1.8',
            re.compile(
                r'(?P<prefix>[\w-]+:\d+:\d+:\d+:\d+#([ACGTN|\d])+)'
                r'('
                r'(?P<linker>/)'
                r'(?P<mate>[12])'
                r'(?P<suffix>)'
                r')?'
            )
        ),
    ],
)


class ReadLength(BaseModel):
    """Read length of a sequencing file.

    Args:
        min: Minimum read length.
        max: Maximum read length.
        mean: Mean of read lengths.
        median: Median of read lengths.
        mode: Mode of read length.

    Attributes:
        min: Minimum read length.
        max: Maximum read length.
        mean: Mean of read lengths.
        median: Median of read lengths.
        mode: Mode of read length.
    """
    min: Optional[int] = None
    max: Optional[int] = None
    mean: Optional[float] = None
    median: Optional[float] = None
    mode: Optional[int] = None


class Stats(BaseModel):
    """Library statistics of an individual sequencing file.

    Args:
        read_length: Tuple of mininimum, maximum, mean, median and mode
            of lengths of reads in library.

    Attributes:
        read_length: Tuple of mininimum, maximum, mean, median and mode
            of lengths of reads in library.
    """
    read_length: ReadLength = ReadLength()


class ResultsStats(BaseModel):
    """Container class for aggregating library statistics information.

    Args:
        file_1: Library statistics for the first file.
        file_2: Library statistics for the second file.

    Attributes:
        file_1: Library statistics for the first file.
        file_2: Library statistics for the second file.
    """
    file_1: Stats = Stats()
    file_2: Stats = Stats()


class StatesType(Enum):
    """Possible outcomes of determining the sequencing library type of an
    individual FASTQ file.

    Attributes:
        file_problem: There was a problem with opening or parsing the file.
        first_mate: All of the sequence identifiers of the processed file
            counts indicate that the library represents the first mate of a
            paired-end library.
        mixed_mates: All of the sequence identifiers of the processed file
            include mate information. However, the file includes at least one
            record for either mate, indicating that the library represents a
            mixed mate library.
        not_available: Library type information is not available for a given
            file, either because no file was provided, the file could not be
            parsed, a library type has not yet been assigned, the processed
            file contains records with sequence identifiers of an unknown
            format, of different formats or that are inconsistent in that they
            indicate the library represents both a single-ended and
            paired-ended library at the same time.
        second_mate: All of the sequence identifiers of the processed file
            indicate that the library represents the second mate of a
            paired-end library.
        single: All of the sequence identifiers of the processed file indicate
            that the library represents a single-end library.
    """
    first_mate = "first_mate"
    mixed_mates = "mixed_mates"
    not_available = None
    second_mate = "second_mate"
    single = "single"


class StatesTypeRelationship(Enum):
    """Possible outcomes of determining the sequencing library type/mate
    relationship between two FASTQ files.

    Attributes:
        not_available: Mate relationship information is not available, likely
            because only a single file was provided or because the mate
            relationship has not yet been evaluated.
        not_mates: The library type information of the files is not compatible,
            either because not a pair of first and second mate files was
            provided, or because the files do not compatible sequence
            identifiers.
        split_mates: One of the provided files represents the first and the
            other the second mates of a paired-end library.
    """
    split_mates = "split_mates"
    not_available = None
    not_mates = "not_mates"


class ResultsType(BaseModel):
    """Container class for aggregating library type and mate relationship
    information.

    Args:
        file_1: Library type of the first file.
        file_2: Library type of the second file.
        relationship: Type/mate relationship between the provided files.

    Attributes:
        file_1: Library type of the first file.
        file_2: Library type of the second file.
        relationship: Type/mate relationship between the provided files.
    """
    file_1: Optional[StatesType] = StatesType.not_available
    file_2: Optional[StatesType] = StatesType.not_available
    relationship: Optional[StatesTypeRelationship] = (
        StatesTypeRelationship.not_available
    )


class Source(BaseModel):
    """Library source of an individual sequencing file.

    Args:
        short_name: Library source short name, e.g., "hsapiens".
        taxon_id: Library source taxon identifer, e.g., `9606`.

    Attributes:
        short_name: Library source short name, e.g., "hsapiens".
        taxon_id: Library source taxon identifer, e.g., `9606`.
    """
    short_name: Optional[str] = None
    taxon_id: Optional[int] = None


class ResultsSource(BaseModel):
    """Container class for aggregating library source.

    Args:
        file_1: Library source of the first file.
        file_2: Library source of the second file.

    Attributes:
        file_1: Library source of the first file.
        file_2: Library source of the second file.
    """
    file_1: Source = Source()
    file_2: Source = Source()


class StatesOrientation(Enum):
    """Enumerator of read orientation types for individual library files. Cf.
    https://salmon.readthedocs.io/en/latest/library_type.html

    Attributes:
        not_available: Orientation type information is not available for a
            given file, either because no file was provided, the file could not
            be parsed, an orientation type has not yet been assigned.
        stranded_forward: Reads are stranded and come from the forward strand.
        stranded_reverse: Reads are stranded and come from the reverse strand.
        unstranded: Reads are unstranded.
    """
    not_available = None
    stranded_forward = "SF"
    stranded_reverse = "SR"
    unstranded = "U"


class StatesOrientationRelationship(Enum):
    """Enumerator of read orientation type relationships for paired-ended
    libraries. Cf. https://salmon.readthedocs.io/en/latest/library_type.html

    Attributes:
        inward_stranded_forward: Mates are oriented toward each other, the
            library is stranded, and first mates come from the forward strand.
        inward_stranded_reverse: Mates are oriented toward each other, the
            library is stranded, and first mates come from the reverse strand.
        inward_unstranded: Mates are oriented toward each other and the library
            is unstranded.
        not_available: Orientation type relationship information is not
            available, likely because only a single file was provided or
            because the orientation type relationship has not been or could not
            be evaluated.
    """
    inward_stranded_forward = "ISF"
    inward_stranded_reverse = "ISR"
    inward_unstranded = "IU"
    not_available = None


class ResultsOrientation(BaseModel):
    """Container class for aggregating library orientation.
     Args:
        file_1: Read orientation of first file.
        file_2: Read orientation of second file.
        relationship: Orientation type relationship between the provided files.

    Attributes:
        file_1: Read orientation of first file.
        file_2: Read orientation of second file.
    """
    file_1: Optional[StatesOrientation] = StatesOrientation.not_available
    file_2: Optional[StatesOrientation] = StatesOrientation.not_available
    relationship: Optional[StatesOrientationRelationship] = (
        StatesOrientationRelationship.not_available
    )


class Layout(BaseModel):
    """Read layout of a single sequencing file.

    Args:
        adapt_3: Adapter sequence ligated to 3'-end of sequence.
        polyA_frac: Fraction of reads containing polyA tails.

    Attributes:
        adapt_3: Adapter sequence ligated to 3'-end of sequence.
        polyA_frac: Fraction of reads containing polyA tails.
    """
    adapt_3: Optional[str] = None
    polyA_frac: Optional[float] = None


class ResultsLayout(BaseModel):
    """Container class for read layout of a sequencing library.

    Args:
        file_1: Adapter sequence present in first file.
        file_2: Adapter sequence present in second file.

    Attributes:
        file_1: Adapter sequence present in first file.
        file_2: Adapter sequence present in second file.
    """
    file_1: Layout = Layout()
    file_2: Layout = Layout()


class Results(BaseModel):
    """Container class for aggregating results from the different inference
    functionalities.

    Args:
        library_type: Library type inference results.
        library_source: Library source inference results.
        orientation: Read orientation inference results.
        read_layout: Read layout inference results.

    Args:
        type: Library type inference results.
        source: Library source inference results.
        read_orientation: Read orientation inference results.
        read_layout: Read layout inference results.
    """
    library_stats: ResultsStats = ResultsStats()
    library_type: ResultsType = ResultsType()
    library_source: ResultsSource = ResultsSource()
    read_orientation: ResultsOrientation = ResultsOrientation()
    read_layout: ResultsLayout = ResultsLayout()


class Args(BaseModel):
    """Configuration model for CLI arguments.

    Args:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        cleanup_regime: Which data to keep after run concludes; one of
            `CleanupRegimes`.
        records: Number of input file records to process; set to `0` to
            process all records.
        threads: Number of threads to run STAR with.
        transcripts_file: File path to transcripts FASTA file.
        read_layout_adapter_file: Path to text file containing 3' adapter
            sequences to scan for (one sequence per line).
        read_layout_min_match_pct: Minimum percentage of reads that contain a
            given adapter in order for it to be considered as the library's
            3'-end adapter.
        read_layout_min_freq_ratio: Minimum frequency ratio between the first
            and second most frequent adapter in order for the former to be
            considered as the library's 3'-end adapter.
        lib_source_min_match_pct: Minimum percentage of reads that are
            consistent with a given source in order for it to be considered as
            the to be considered the library's source.
        lib_source_min_freq_ratio: Minimum frequency ratio between the first
            and second most frequent source in order for the former to be
            considered the library's source.
        read_orientation_min_mapped_reads: Minimum number of mapped reads for
            deeming the read orientation result reliable.
        read_orientation_min_fraction: Minimum fraction of mapped reads
            required to be consistent with a given read orientation state in
            order for that orientation to be reported. Must be above 0.5.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        out_dir: Path to directory where output is written to.
        run_id: Random string identifier for HTSinfer run.
        tmp_dir: Path to directory where temporary output is written to.
        cleanup_regime: Which data to keep after run concludes; one of
            `CleanupRegimes`.
        records: Number of input file records to process.
        threads: Number of threads to run STAR with.
        transcripts_file: File path to transcripts FASTA file.
        read_layout_adapter_file: Path to text file containing 3' adapter
            sequences to scan for (one sequence per line).
        read_layout_min_match_pct: Minimum percentage of reads that contain a
            given adapter in order for it to be considered as the library's
            3'-end adapter.
        read_layout_min_freq_ratio: Minimum frequency ratio between the first
            and second most frequent adapter in order for the former to be
            considered as the library's 3'-end adapter.
        lib_source_min_match_pct: Minimum percentage of reads that are
            consistent with a given source in order for it to be considered as
            the to be considered the library's source.
        lib_source_min_freq_ratio: Minimum frequency ratio between the first
            and second most frequent source in order for the former to be
            considered the library's source.
        lib_type_max_distance: Upper limit on the difference in the
            reference sequence coordinates between the two mates to be
            considered as coming from a single fragment.
            (Used only when sequence identifiers do not match)
        lib_type_mates_cutoff: Minimum fraction of mates that can be mapped to
            compatible loci and are considered concordant pairs / all mates.
            (Used only when sequence identifiers do not match)
        read_orientation_min_mapped_reads: Minimum number of mapped reads for
            deeming the read orientation result reliable.
        read_orientation_min_fraction: Minimum fraction of mapped reads
            required to be consistent with a given read orientation state in
            order for that orientation to be reported. Must be above 0.5.
        path_1_processed: Path to processed `path_1` file.
        path_2_processed: Path to processed `path_2` file.
        t_file_processed: Path to processed `transcripts_file` file.
        state: State of the run; one of `RunStates`.
        results: Results container for storing determined library metadata.
    """
    path_1: Path = Path()
    path_2: Optional[Path] = None
    out_dir: Path = Path.cwd() / 'results_htsinfer'
    tmp_dir: Path = \
        Path(tempfile.gettempdir()) / 'tmp_htsinfer'
    cleanup_regime: CleanupRegimes = \
        CleanupRegimes.DEFAULT
    records: int = 1000000
    threads: int = 1
    transcripts_file: Path = Path()
    read_layout_adapter_file: Path = Path()
    read_layout_min_match_pct: float = 0.1
    read_layout_min_freq_ratio: float = 2
    lib_source_min_match_pct: float = 2
    lib_source_min_freq_ratio: float = 2
    lib_type_max_distance: int = 1000
    lib_type_mates_cutoff: float = 0.95
    read_orientation_min_mapped_reads: int = 20
    read_orientation_min_fraction: float = 0.75
    path_1_processed: Path = Path()
    path_2_processed: Optional[Path] = None
    t_file_processed: Path = Path()


class Config(BaseModel):
    """Configuration model for CLI arguments and inference results.

    Args:
        args: Container class for CLI arguments.
        results: Container class for aggregating results
                from the different inference functionalities.
    Attributes:
        args: Container class for CLI arguments.
        results: Container class for aggregating results
                from the different inference functionalities.
    """
    args: Args = Args()
    results: Results = Results()
