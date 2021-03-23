"""Main module."""
# pylint: disable=fixme

from enum import (Enum, IntEnum)
import logging
from pathlib import Path
from random import choices
import shutil
import string
import sys
import tempfile
from typing import Optional

from htsinfer.exceptions import WorkEnvProblem
from htsinfer.models import Results

LOGGER = logging.getLogger(__name__)


class CleanupRegimes(Enum):
    """Enumerator of cleanup regimes."""
    default = "default"
    keep_all = "keep_all"
    keep_none = "keep_none"
    keep_results = "keep_results"


class RunStates(IntEnum):
    """Enumerator of run states and exit codes."""
    okay = 0
    warning = 1
    error = 2


class HtsInfer:  # pylint: disable=too-many-instance-attributes
    """Determine sequencing library metadata.

    Args:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        cleanup_regime: Which data to keep after run concludes; one of

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        out_dir: Path to directory where output is written to.
        run_id: Random string identifier for HTSinfer run.
        tmp_dir: Path to directory where temporary output is written to.
        cleanup_regime: Which data to keep after run concludes; one of
            `CleanupRegimes`.
        path_1_processed: Path to processed `path_1` file.
        path_2_processed: Path to processed `path_2` file.
        state: State of the run; one of `RunStates`.
        results: Results container for storing determined library metadata.
    """

    def __init__(  # pylint: disable=too-many-arguments
        self,
        path_1: Path,
        path_2: Optional[Path] = None,
        out_dir: Path = Path.cwd(),
        tmp_dir: Path = Path(tempfile.gettempdir()),
        cleanup_regime: CleanupRegimes = CleanupRegimes.default,
    ):
        """Class constructor."""
        self.path_1 = path_1
        self.path_2 = path_2
        self.run_id = ''.join(
            choices(string.ascii_uppercase + string.digits, k=5)
        )
        self.out_dir = out_dir / self.run_id
        self.tmp_dir = tmp_dir / f"tmp_{self.run_id}"
        self.cleanup_regime = cleanup_regime
        self.path_1_processed: Path = self.path_1
        self.path_2_processed: Optional[Path] = self.path_2
        self.state: RunStates = RunStates.okay
        self.results: Results = Results()

    def evaluate(self):
        """Determine library metadata."""
        try:
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

        except WorkEnvProblem as exc:
            self.state = RunStates.error
            LOGGER.error(f"{type(exc).__name__}: {str(exc)}")

        # log results
        LOGGER.info(f"Results: {self.results.json()}")

    def prepare_env(self):
        """Set up work environment."""
        # create results directory
        try:
            self.out_dir.mkdir()
        except OSError:
            raise WorkEnvProblem(
                f"Creation of results directory failed: {self.out_dir}"
            )
        LOGGER.info(f"Created results directory: {self.out_dir}")
        # create temporary directory
        try:
            self.tmp_dir.mkdir()
        except OSError:
            raise WorkEnvProblem(
                f"Creation of temporary directory failed: {self.tmp_dir}"
            )
        LOGGER.info(f"Created temporary directory: {self.tmp_dir}")
        LOGGER.debug("Created work environment")

    def process_inputs(self):
        """Process and validate inputs."""
        # TODO: implement

    def get_library_type(self):
        """Determine library type."""
        # TODO: implement

    def get_library_source(self):
        """Determine library source."""
        # TODO: implement

    def get_read_orientation(self):
        """Determine read orientation."""
        # TODO: implement

    def get_read_layout(self):
        """Determine read layout."""
        # TODO: implement

    def clean_up(self):
        """Clean up work environment."""
        # set default cleanup regime
        if self.cleanup_regime is CleanupRegimes.default:
            if (
                logging.root.level == logging.DEBUG or
                self.state is RunStates.error
            ):
                self.cleanup_regime = CleanupRegimes.keep_all
            elif self.state is RunStates.warning:
                self.cleanup_regime = CleanupRegimes.keep_results
            else:
                self.cleanup_regime = CleanupRegimes.keep_none
        LOGGER.debug(f"Cleanup regime: {self.cleanup_regime}")

        # remove results directory
        if self.cleanup_regime == CleanupRegimes.keep_none:
            try:
                shutil.rmtree(self.out_dir)
            except OSError:
                raise WorkEnvProblem(
                    f"Removal of results directory failed: {self.out_dir}"
                )
            LOGGER.info(f"Removed results directory: {self.out_dir}")

        # remove temporary directory
        if (
            self.cleanup_regime == CleanupRegimes.keep_results or
            self.cleanup_regime == CleanupRegimes.keep_none
        ):
            try:
                shutil.rmtree(self.tmp_dir)
            except OSError:
                raise WorkEnvProblem(
                    f"Removal of temporary directory failed: {self.tmp_dir}"
                )
            LOGGER.info(f"Removed temporary directory: {self.tmp_dir}")

    def print(self):
        """Print results to STDOUT."""
        sys.stdout.write(self.results.json())
