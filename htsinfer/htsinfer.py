"""Main module."""

from enum import IntEnum
import logging
from pathlib import Path
import sys
from typing import Optional

from htsinfer.models import Results

LOGGER = logging.getLogger(__name__)


class RunStates(IntEnum):
    """Enumerator of run states and exit codes."""
    okay = 0
    warning = 1
    error = 2


class HtsInfer:
    """Determine sequencing library metadata.

    Args:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        state: State of the run; one of `RunStates`.
        results: Results container for storing determined library metadata.
    """

    def __init__(
        self,
        path_1: Path,
        path_2: Optional[Path] = None,
    ):
        """Class constructor."""
        self.path_1 = path_1
        self.path_2 = path_2
        self.state: RunStates = RunStates.okay
        self.results: Results = Results()

    def evaluate(self):
        """Determine library metadata."""
        # set up work environment
        LOGGER.info("Setting up work environment...")
        self.prepare_env()

        # preprocess inputs
        LOGGER.info("Processing and validating input data...")
        self.process_inputs()

        # determine library type
        LOGGER.info("Determining library type...")
        self.get_library_type()

        # determine library source
        LOGGER.info("Determining library source...")
        self.get_library_source()

        # determine read orientation
        LOGGER.info("Determining read orientation...")
        self.get_read_orientation()

        # determine read layout
        LOGGER.info("Determining read layout...")
        self.get_read_layout()

        # postprocessing
        LOGGER.info("Cleaning up work environment...")
        self.clean_up()

        # log results
        LOGGER.info(f"Results: {self.results.json()}")

    def prepare_env(self):
        """Prepare work environment."""
        # TODO: implement  # pylint: disable=fixme

    def process_inputs(self):
        """Process and validate inputs."""
        # TODO: implement  # pylint: disable=fixme

    def get_library_type(self):
        """Determine library type."""
        # TODO: implement  # pylint: disable=fixme

    def get_library_source(self):
        """Determine library source."""
        # TODO: implement  # pylint: disable=fixme

    def get_read_orientation(self):
        """Determine read orientation."""
        # TODO: implement  # pylint: disable=fixme

    def get_read_layout(self):
        """Determine read layout."""
        # TODO: implement  # pylint: disable=fixme

    def clean_up(self):
        """Clean up working environment."""
        # TODO: implement  # pylint: disable=fixme

    def print(self):
        """Print results to STDOUT."""
        sys.stdout.write(self.results.json())
