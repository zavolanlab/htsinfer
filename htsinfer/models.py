"""Data models."""
# pylint: disable=too-few-public-methods

from typing import Optional

from pydantic import BaseModel  # pylint: disable=no-name-in-module


class ResultsType(BaseModel):
    """TODO: implement"""


class ResultsSource(BaseModel):
    """TODO: implement"""


class ResultsReadOrientation(BaseModel):
    """TODO: implement"""


class ResultsReadLayout(BaseModel):
    """TODO: implement"""


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
    library_type: Optional[ResultsType] = None
    library_source: Optional[ResultsSource] = None
    read_orientation: Optional[ResultsReadOrientation] = None
    read_layout: Optional[ResultsReadLayout] = None
