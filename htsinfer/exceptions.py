"""Custom exceptions."""


class FileProblem(Exception):
    """Exception raised when file could not be opened or parsed."""


class WorkEnvProblem(Exception):
    """Exception raised when the work environment could not be set up or
    cleaned."""
