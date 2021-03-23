"""Custom exceptions."""


class FileProblem(Exception):
    """Exception raised when file could not be opened or parsed."""


class InconsistentFastqIdentifiers(Exception):
    """Exception raised when inconsistent FASTQ sequence identifiers were
    ecountered.
    """


class MetadataWarning(Exception):
    """Exception raised when metadata could not be determined."""


class UnknownFastqIdentifier(Exception):
    """Exception raised when a FASTQ sequence identifier of unknown format was
    ecountered.
    """


class WorkEnvProblem(Exception):
    """Exception raised when the work environment could not be set up or
    cleaned."""
