"""Main module."""
# pylint: disable=fixme

import logging
from pathlib import Path
from random import choices
import shutil
import string
import sys
import tempfile
from typing import Optional

from htsinfer.exceptions import (
    FileProblem,
    MetadataWarning,
    WorkEnvProblem,
)
from htsinfer.get_library_type import GetLibType
from htsinfer.models import (
    CleanupRegimes,
    Results,
    RunStates,
)
from htsinfer.subset_fastq import SubsetFastq

LOGGER = logging.getLogger(__name__)


class HtsInfer:  # pylint: disable=too-many-instance-attributes
    """Determine sequencing library metadata.

    Args:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        out_dir: Path to directory where output is written to.
        tmp_dir: Path to directory where temporary output is written to.
        cleanup_regime: Which data to keep after run concludes; one of
        records: Number of input file records to process; set to `0` to
            process all records.
            `CleanupRegimes`.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        out_dir: Path to directory where output is written to.
        run_id: Random string identifier for HTSinfer run.
        tmp_dir: Path to directory where temporary output is written to.
        cleanup_regime: Which data to keep after run concludes; one of
            `CleanupRegimes`.
        records: Number of input file records to process.
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
        cleanup_regime: CleanupRegimes = CleanupRegimes.DEFAULT,
        records: int = 0,
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
        self.records = records
        self.path_1_processed: Path = self.path_1
        self.path_2_processed: Optional[Path] = self.path_2
        self.state: RunStates = RunStates.OKAY
        self.results: Results = Results()

    def evaluate(self):
        """Determine library metadata."""
        try:
            # set up work environment
            LOGGER.info("Setting up work environment...")
            self.prepare_env()

            try:
                # preprocess inputs
                LOGGER.info("Processing and validating input data...")
                self.process_inputs()

                # determine library type
                LOGGER.info("Determining library type...")
                try:
                    self.get_library_type()
                except MetadataWarning as exc:
                    self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")
                LOGGER.info(
                    "Library type determined: "
                    f"{self.results.library_type.json()}"
                )

                # determine library source
                LOGGER.info("Determining library source...")
                try:
                    self.get_library_source()
                except MetadataWarning as exc:
                    self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")

                # determine read orientation
                LOGGER.info("Determining read orientation...")
                try:
                    self.get_read_orientation()
                except MetadataWarning as exc:
                    self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")

                # determine read layout
                LOGGER.info("Determining read layout...")
                try:
                    self.get_read_layout()
                except MetadataWarning as exc:
                    self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")

            except FileProblem as exc:
                self.state = RunStates.ERROR
                LOGGER.error(f"{type(exc).__name__}: {str(exc)}")

            # postprocessing
            LOGGER.info("Cleaning up work environment...")
            self.clean_up()

        except WorkEnvProblem as exc:
            self.state = RunStates.ERROR
            LOGGER.error(f"{type(exc).__name__}: {str(exc)}")

        # log results
        LOGGER.info(f"Results: {self.results.json()}")

    def prepare_env(self):
        """Set up work environment."""
        # create results directory
        try:
            self.out_dir.mkdir()
        except OSError as exc:
            raise WorkEnvProblem(
                f"Creation of results directory failed: {self.out_dir}"
            ) from exc
        LOGGER.info(f"Created results directory: {self.out_dir}")

        # create temporary directory
        try:
            self.tmp_dir.mkdir()
        except OSError as exc:
            raise WorkEnvProblem(
                f"Creation of temporary directory failed: {self.tmp_dir}"
            ) from exc
        LOGGER.info(f"Created temporary directory: {self.tmp_dir}")
        LOGGER.debug("Created work environment")

    def process_inputs(self):
        """Process and validate inputs."""
        # process first file
        input_files_1 = SubsetFastq(
            path=self.path_1,
            out_dir=self.tmp_dir,
            records=self.records,
        )
        input_files_1.process()
        self.path_1_processed = input_files_1.out_path
        LOGGER.info(f"Location processed file 1: {self.path_1_processed}")

        # process second file, if available
        if self.path_2 is not None:
            input_files_2 = SubsetFastq(
                path=self.path_2,
                out_dir=self.tmp_dir,
                records=self.records,
            )
            input_files_2.process()
            self.path_2_processed = input_files_2.out_path

    def get_library_type(self):
        """Determine library type."""
        get_lib_type = GetLibType(
            path_1=self.path_1_processed,
            path_2=self.path_2_processed,
        )
        get_lib_type.evaluate()
        self.results.library_type = get_lib_type.results

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
        if self.cleanup_regime is CleanupRegimes.DEFAULT:
            if (
                logging.root.level == logging.DEBUG or
                self.state is RunStates.ERROR
            ):
                self.cleanup_regime = CleanupRegimes.KEEP_ALL
            elif self.state is RunStates.WARNING:
                self.cleanup_regime = CleanupRegimes.KEEP_RESULTS
            else:
                self.cleanup_regime = CleanupRegimes.KEEP_NONE
        LOGGER.debug(f"Cleanup regime: {self.cleanup_regime}")

        # remove results directory
        if self.cleanup_regime == CleanupRegimes.KEEP_NONE:
            try:
                shutil.rmtree(self.out_dir)
            except OSError as exc:
                raise WorkEnvProblem(
                    f"Removal of results directory failed: {self.out_dir}"
                ) from exc
            LOGGER.info(f"Removed results directory: {self.out_dir}")

        # remove temporary directory
        if (
            self.cleanup_regime == CleanupRegimes.KEEP_RESULTS or
            self.cleanup_regime == CleanupRegimes.KEEP_NONE
        ):
            try:
                shutil.rmtree(self.tmp_dir)
            except OSError as exc:
                raise WorkEnvProblem(
                    f"Removal of temporary directory failed: {self.tmp_dir}"
                ) from exc
            LOGGER.info(f"Removed temporary directory: {self.tmp_dir}")

    def print(self):
        """Print results to STDOUT."""
        sys.stdout.write(self.results.json())
