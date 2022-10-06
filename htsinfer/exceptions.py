"""Custom exceptions."""


class FileProblem(Exception):
    """Exception raised when file could not be opened or parsed."""


class InconsistentFastqIdentifiers(Exception):
    """Exception raised when inconsistent FASTQ sequence identifiers were
    ecountered.
    """


class KallistoProblem(Exception):
    """Exception raised when running kallisto index and quant commands."""


class MetadataWarning(Exception):
    """Exception raised when metadata could not be determined."""


class StarProblem(Exception):
    """Exception raised when running STAR index and quant commands."""


class UnknownFastqIdentifier(Exception):
    """Exception raised when a FASTQ sequence identifier of unknown format was
    ecountered.
    """


class WorkEnvProblem(Exception):
    """Exception raised when the work environment could not be set up or
    cleaned."""


class TranscriptsFastaProblem(Exception):
    """Exception raised when an invalid transcripts fasta file is passed."""
