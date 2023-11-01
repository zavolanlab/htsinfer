"""Main module."""

from functools import partial
import gzip
import logging
from os import linesep
from pathlib import Path
from random import choices
import shutil
import string
import sys

from htsinfer.exceptions import (
    FileProblem,
    MetadataWarning,
    WorkEnvProblem,
)
from htsinfer.get_library_source import GetLibSource
from htsinfer.get_library_stats import GetLibStats
from htsinfer.get_library_type import GetLibType
from htsinfer.get_read_orientation import GetOrientation
from htsinfer.get_read_layout import GetReadLayout
from htsinfer.models import (
    CleanupRegimes,
    ResultsSource,
    RunStates,
    Config,
)
from htsinfer.subset_fastq import SubsetFastq
from htsinfer.mapping import Mapping

LOGGER = logging.getLogger(__name__)


class HtsInfer:
    """Determine sequencing library metadata.

    Args:
        config: Container class for all arguments used in inference
                and results produced by the class.

    Attributes:
        config: Container class for all arguments used in inference
                and results produced by the class.
        run_id: Random string identifier for HTSinfer run.
        state: State of the run; one of `RunStates`.
    """
    def __init__(
        self,
        config: Config,
    ):
        """Class constructor."""
        self.config = config
        self.run_id = ''.join(
            choices(string.ascii_uppercase + string.digits, k=5)
        )
        self.config.args.out_dir = Path(config.args.out_dir) / self.run_id
        self.config.args.tmp_dir = \
            Path(config.args.tmp_dir) / f"tmp_{self.run_id}"
        self.config.args.path_1_processed = self.config.args.path_1
        self.config.args.path_2_processed = self.config.args.path_2
        self.config.args.t_file_processed = (
            config.args.tmp_dir / config.args.transcripts_file.stem
            if config.args.transcripts_file.suffix == ".gz"
            else config.args.tmp_dir / config.args.transcripts_file.name
        )
        self.state: RunStates = RunStates.OKAY
        self.mapping: Mapping = Mapping(config=self.config)

    def evaluate(self):
        """Determine library metadata."""
        try:
            # set up work environment
            LOGGER.info("Setting up work environment...")
            LOGGER.info(f"Run identifier: {self.run_id}")
            self.prepare_env()

            try:
                # preprocess inputs
                LOGGER.info("Processing and validating input data...")
                self.process_inputs()

                # determine library stats
                LOGGER.info("Determining library statistics...")
                self.get_library_stats()
                LOGGER.info(
                    "Library stats determined: "
                    f"{self.config.results.library_stats.json()}"
                )

                # determine library source
                LOGGER.info("Determining library source...")
                self.config.results.library_source = self.get_library_source()
                LOGGER.info(
                    "Library source determined: "
                    f"{self.config.results.library_source.json()}"
                )

                # determine library type
                LOGGER.info("Determining library type...")
                try:
                    self.get_library_type()
                except MetadataWarning as exc:
                    if self.state is RunStates.OKAY:
                        self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")
                LOGGER.info(
                    "Library type determined: "
                    f"{self.config.results.library_type.json()}"
                )

                # determine read orientation
                LOGGER.info("Determining read orientation...")
                try:
                    self.get_read_orientation()
                except MetadataWarning as exc:
                    if self.state is RunStates.OKAY:
                        self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")
                LOGGER.info(
                    "Read orientation determined: "
                    f"{self.config.results.read_orientation.json()}"
                )

                # determine read layout
                LOGGER.info("Determining read layout...")
                try:
                    self.get_read_layout()
                except MetadataWarning as exc:
                    if self.state is RunStates.OKAY:
                        self.state = RunStates.WARNING
                    LOGGER.warning(f"{type(exc).__name__}: {str(exc)}")
                LOGGER.info(
                    "Read layout determined: "
                    f"{self.config.results.read_layout.json()}"
                )

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
        LOGGER.info(f"Results: {self.config.results.json()}")

    def prepare_env(self):
        """Set up work environment."""
        # create results directory
        try:
            self.config.args.out_dir.mkdir(parents=True)
        except OSError as exc:
            raise WorkEnvProblem(
                f"Creation of results directory failed: "
                f"{self.config.args.out_dir}"
            ) from exc
        LOGGER.info(f"Created results directory: {self.config.args.out_dir}")

        # create temporary directory
        try:
            self.config.args.tmp_dir.mkdir(parents=True)
        except OSError as exc:
            raise WorkEnvProblem(
                f"Creation of temporary directory failed: "
                f"{self.config.args.tmp_dir}"
            ) from exc
        LOGGER.info(f"Created temporary directory: {self.config.args.tmp_dir}")
        LOGGER.debug("Created work environment")

    def process_inputs(self):
        """Process and validate inputs."""

        # validate input parameters
        if self.config.args.read_orientation_min_fraction <= 0.5:
            raise ValueError(
                "Value for parameter 'read_orientation_min_fraction' outside"
                "permitted boundaries; expected: >0.5; found: "
                f"{self.config.args.read_orientation_min_fraction}"
            )

        # process first file
        LOGGER.debug(f"Processing read file 1: {self.config.args.path_1}")
        input_files_1 = SubsetFastq(
            path=self.config.args.path_1,
            out_dir=self.config.args.tmp_dir,
            records=self.config.args.records,
        )
        input_files_1.process()
        self.config.args.path_1_processed = input_files_1.out_path
        LOGGER.info(f"Processed read file 1: "
                    f"{self.config.args.path_1_processed}")

        # process second file, if available
        if self.config.args.path_2 is not None:
            LOGGER.debug(f"Processing read file 2: {self.config.args.path_2}")
            input_files_2 = SubsetFastq(
                path=self.config.args.path_2,
                out_dir=self.config.args.tmp_dir,
                records=self.config.args.records,
            )
            input_files_2.process()
            self.config.args.path_2_processed = input_files_2.out_path
            LOGGER.info(f"Processed read file 2: "
                        f"{self.config.args.path_2_processed}")

        # process transcripts file
        LOGGER.debug(f"Processing transcripts file: "
                     f"{self.config.args.transcripts_file}")
        if self.config.args.transcripts_file.suffix == ".gz":
            _open = partial(gzip.open)
        else:
            _open = open
        try:
            with _open(self.config.args.transcripts_file, 'rb') as f_in:
                with open(self.config.args.t_file_processed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        except Exception as exc:
            raise FileProblem(exc) from exc
        LOGGER.info(
            f"Processed transcripts file: {self.config.args.t_file_processed}"
        )

    def get_library_stats(self):
        """Determine library statistics."""
        get_lib_stats = GetLibStats(
            config=self.config,
        )
        self.config.results.library_stats = get_lib_stats.evaluate()

    def get_library_source(self) -> ResultsSource:
        """Determine library source.

        Returns:
            Library source results.
        """
        get_lib_source = GetLibSource(
            config=self.config,
        )
        results = get_lib_source.evaluate()
        return results

    def get_library_type(self):
        """Determine library type."""
        get_lib_type = GetLibType(
            config=self.config,
            mapping=self.mapping,
        )
        get_lib_type.evaluate()
        self.config.results.library_type = get_lib_type.results

    def get_read_orientation(self):
        """Determine read orientation."""
        get_read_orientation = GetOrientation(
            config=self.config,
            mapping=self.mapping,
        )
        self.config.results.read_orientation = get_read_orientation.evaluate()

    def get_read_layout(self):
        """Determine read layout."""
        get_read_layout = GetReadLayout(
            config=self.config,
        )
        get_read_layout.evaluate()
        self.config.results.read_layout = get_read_layout.results

    def clean_up(self):
        """Clean up work environment."""
        # set default cleanup regime
        if self.config.args.cleanup_regime is CleanupRegimes.DEFAULT:
            if (
                logging.root.level == logging.DEBUG or
                self.state is RunStates.ERROR
            ):
                self.config.args.cleanup_regime = CleanupRegimes.KEEP_ALL
            elif self.state is RunStates.WARNING:
                self.config.args.cleanup_regime = CleanupRegimes.KEEP_RESULTS
            else:
                self.config.args.cleanup_regime = CleanupRegimes.KEEP_NONE
        LOGGER.debug(f"Cleanup regime: {self.config.args.cleanup_regime}")

        # remove results directory
        if self.config.args.cleanup_regime == CleanupRegimes.KEEP_NONE:
            try:
                shutil.rmtree(self.config.args.out_dir)
            except OSError as exc:
                raise WorkEnvProblem(
                    f"Removal of results directory failed: "
                    f"{self.config.args.out_dir}"
                ) from exc
            LOGGER.info(f"Removed results directory: "
                        f"{self.config.args.out_dir}")

        # remove temporary directory
        if self.config.args.cleanup_regime in (
            CleanupRegimes.KEEP_RESULTS,
            CleanupRegimes.KEEP_NONE
        ):
            try:
                shutil.rmtree(self.config.args.tmp_dir)
            except OSError as exc:
                raise WorkEnvProblem(
                    f"Removal of temporary directory failed: "
                    f"{self.config.args.tmp_dir}"
                ) from exc
            LOGGER.info(f"Removed temporary directory: "
                        f"{self.config.args.tmp_dir}")

    def print(self):
        """Print results to STDOUT."""
        sys.stdout.write(
            self.config.results.model_dump_json(
                indent=3,
            ) + linesep
        )
