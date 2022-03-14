"""Data models."""

from enum import (
    Enum,
    IntEnum,
)
import logging
import re
from typing import Optional

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
        matching_stranded_forward: Mates are oriented in the same direction,
            the library is stranded, and both mates come from the forward
            strand.
        matching_stranded_reverse: Mates are oriented in the same direction,
            the library is stranded, and both mates come from the reverse
            strand.
        matching_unstranded: Mates are oriented in the same direction and the
            library is unstranded.
        not_available: Orientation type relationship information is not
            available, likely because only a single file was provided or
            because the orientation type relationship has not yet been
            evaluated.
        outward_stranded_forward: Mates are oriented away from each other, the
            library is stranded, and first mates come from the forward strand.
        outward_stranded_reverse: Mates are oriented away from each other, the
            library is stranded, and first mates come from the reverse strand.
        outward_unstranded: Mates are oriented away from each other and the
            library is unstranded.
    """
    inward_stranded_forward = "ISF"
    inward_stranded_reverse = "ISR"
    inward_unstranded = "IU"
    matching_stranded_forward = "MSF"
    matching_stranded_reverse = "MSR"
    matching_unstranded = "MU"
    not_available = None
    outward_stranded_forward = "OSF"
    outward_stranded_reverse = "OSR"
    outward_unstranded = "OU"


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


class StatesType(Enum):
    """Possible outcomes of determining the sequencing library type of an
    individual FASTQ file.

    Attributes:
        file_problem: There was a problem with opening or parsing the file.
        first_mate: All of the sequence identifiers of the processed file
            indicate that the library represents the first mate of a paired-end
            library.
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
            the other the second mates of a paired-end library.
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
    file_1: StatesType = StatesType.not_available
    file_2: StatesType = StatesType.not_available
    relationship: StatesTypeRelationship = (
        StatesTypeRelationship.not_available
    )


class ResultsSource(BaseModel):
    """TODO: implement"""


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
    file_1: StatesOrientation = StatesOrientation.not_available
    file_2: StatesOrientation = StatesOrientation.not_available
    relationship: StatesOrientationRelationship = (
        StatesOrientationRelationship.not_available
    )


class Layout(BaseModel):
    """Read layout of a single sequencing file.

    Args:
        adapt_3: Adapter sequence ligated to 3'-end of sequence.

    Attributes:
        adapt_3: Adapter sequence ligated to 3'-end of sequence.
    """
    adapt_3: Optional[str] = None


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
    library_type: ResultsType = ResultsType()
    library_source: ResultsSource = ResultsSource()
    read_orientation: ResultsOrientation = ResultsOrientation()
    read_layout: ResultsLayout = ResultsLayout()
