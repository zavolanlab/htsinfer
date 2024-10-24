"""Infer read orientation from sample data."""

from collections import defaultdict
import logging
from pathlib import Path
from typing import (Any, DefaultDict, Dict, List)

import pysam  # type: ignore
import pandas as pd  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
)
from htsinfer.models import (
    ResultsOrientation,
    StatesOrientation,
    StatesOrientationRelationship,
    Config,
)
from htsinfer.mapping import Mapping

LOGGER = logging.getLogger(__name__)


class GetOrientation:
    """Determine library strandedness and relative read orientation of a
    single- or paired-end seguencing library.

    Args:
        config: Container class for all arguments used in inference
                and results produced by the class.

    Attributes:
        paths: Tuple of one or two paths for single-end and paired end library
            files.
        library_type: ResultsType object with library type and mate
            relationship.
        library_source: ResultsSource object with source information on each
            library file.
        transcripts_file: File path to an uncompressed transcripts file in
            FASTA format.
        tmp_dir: Path to directory where temporary output is written to.
        threads_star: Number of threads to run STAR with.
        min_mapped_reads: Minimum number of mapped reads for deeming the
            read orientation result reliable.
        min_fraction: Minimum fraction of mapped reads required to be
            consistent with a given read orientation state in order for that
            orientation to be reported. Must be above 0.5.
    """
    def __init__(
        self,
        config: Config,
        mapping: Mapping,
    ):
        """Class contructor."""
        self.paths = (config.args.path_1_processed,
                      config.args.path_2_processed)
        self.library_type = config.results.library_type
        self.library_source = config.results.library_source
        self.transcripts_file = config.args.t_file_processed
        self.tmp_dir = config.args.tmp_dir
        self.out_dir = config.args.out_dir
        self.min_mapped_reads = config.args.read_orientation_min_mapped_reads
        self.min_fraction = config.args.read_orientation_min_fraction
        self.mapping = mapping

    def evaluate(self) -> ResultsOrientation:
        """Infer read orientation.

        Returns:
            Orientation results object.
        """

        self.mapping.paths = self.paths
        self.mapping.library_type = self.library_type
        self.mapping.library_source = self.library_source
        self.mapping.transcripts_file = self.transcripts_file
        self.mapping.tmp_dir = self.tmp_dir

        if not self.mapping.mapped and (
            self.library_source.file_1.short_name is not None or
            self.library_source.file_2.short_name is not None
        ):
            LOGGER.debug("Determining read relationship by alignment...")
            self.mapping.evaluate()

        return self.process_alignments(star_dirs=self.mapping.star_dirs)

    def process_alignments(
        self,
        star_dirs: List[Path],
    ) -> ResultsOrientation:
        """Determine read orientation of one or two single-ended or one
        paired-end sequencing library.

        Args:
            star_dirs: List of one or two paths to STAR output directories.

        Returns:
            Read orientation state of library or libraries.
        """
        results: ResultsOrientation = ResultsOrientation()
        paths = [path / 'Aligned.out.sam' for path in star_dirs]
        if len(star_dirs) == 1:
            if star_dirs[0].name == "paired":
                results = self.process_paired(sam=paths[0])
            if star_dirs[0].name == "file_1":
                results.file_1 = self.process_single(sam=paths[0])
        elif len(star_dirs) == 2:
            results.file_1 = self.process_single(sam=paths[0])
            results.file_2 = self.process_single(sam=paths[1])
        return results

    def process_single(
        self,
        sam: Path,
    ) -> StatesOrientation:
        """Determine read orientation of a single-ended sequencing library.
        Args:
            sam: Path to SAM file.
        Returns:
            Read orientation state of library.

        Raises:
            Sam file could not be processed.
        """
        LOGGER.debug(f"Processing SAM file: '{sam}'")

        # parsing alignments
        states: DefaultDict[str, List[StatesOrientation]] = defaultdict(
            lambda: []
        )
        try:
            with pysam.AlignmentFile(str(sam), "r") as _file:
                for record in _file.fetch():

                    # ensure that "unmapped" flag is not set and query name
                    # is set
                    if (
                        not record.flag & (1 << 2) and
                        isinstance(record.query_name, str)
                    ):

                        # check if reverse complemented
                        if record.flag & (1 << 4):
                            states[record.query_name].append(
                                StatesOrientation.stranded_reverse
                            )
                        elif not record.flag & (1 << 4):
                            states[record.query_name].append(
                                StatesOrientation.stranded_forward
                            )
        except (OSError, ValueError) as exc:
            raise FileProblem(
                f"Failed to open SAM file: '{sam}'"
            ) from exc

        LOGGER.debug("Deciding read orientation...")
        reads = len(states)
        fractions = [
            self.get_frequencies(*state) for state in states.values()
        ]
        fractions_all_states = {
            key: val / reads for (key, val) in
            self.sum_dicts(*fractions).items()
        }
        orientation: StatesOrientation = StatesOrientation.not_available
        if len(fractions_all_states) == 0:
            return orientation
        most_common_state = max(
            fractions_all_states,
            key=fractions_all_states.get  # type: ignore
        )
        fraction_most_common_state = max(fractions_all_states.values())
        if reads >= self.min_mapped_reads:
            if fraction_most_common_state > self.min_fraction:
                orientation = most_common_state
            else:
                orientation = StatesOrientation.unstranded

        orient_df = self.create_orient_df(
            reads, fractions_all_states, orientation, paired=False
        )

        LOGGER.debug(
            f"Required number of mapped reads: {self.min_mapped_reads}"
        )
        LOGGER.debug(f"Number of mapped reads: {orient_df.iloc[0, 0]}")
        LOGGER.debug(f"Fraction of SF: {orient_df.iloc[0, 1]}")
        LOGGER.debug(f"Fraction of SR: {orient_df.iloc[0, 2]}")
        LOGGER.debug(f"Orientation: {orient_df.iloc[0, 3]}")

        self.write_orientation_to_json(orient_df, self.paths[0].name)

        return orientation

    def process_paired(  # pylint: disable=R0912,R0915
        self,
        sam: Path,
    ) -> ResultsOrientation:
        """Determine read orientation of a paired-ended sequencing library.

        Args:
            sam: Path to SAM file.

        Returns:
            Read orientation state of each mate and orientation state
                relationship of library.
        """
        LOGGER.debug(f"Processing SAM file: '{sam}'")

        # parsing alignments
        states: DefaultDict[
            str, List[StatesOrientationRelationship]
        ] = defaultdict(lambda: [])
        try:
            with pysam.AlignmentFile(str(sam), "r") as _file:
                for record_1 in _file.fetch():

                    # ensure both mates are mapped and part of mate pair
                    if not (
                        record_1.flag & (1 << 0) and
                        record_1.flag & (1 << 1) and
                        not record_1.flag & (1 << 2) and
                        not record_1.flag & (1 << 3)
                    ):
                        continue

                    # check which alignment is first mate
                    record_2 = next(_file)
                    if (
                        record_1.flag & (1 << 6) and
                        record_2.flag & (1 << 7)
                    ):
                        mate_1 = record_1
                        mate_2 = record_2
                    elif (
                        record_2.flag & (1 << 6) and
                        record_1.flag & (1 << 7)
                    ):
                        mate_1 = record_2
                        mate_2 = record_1
                    else:
                        continue

                    # check orientation: forward / inward
                    if (
                        not mate_1.flag & (1 << 4) and
                        mate_1.pos < mate_2.pos  # type: ignore
                    ):
                        states[str(mate_1.query_name)].append(
                            StatesOrientationRelationship.
                            inward_stranded_forward
                        )

                    # check orientation: reverse / inward
                    elif (
                        mate_1.flag & (1 << 4) and
                        mate_1.pos > mate_2.pos  # type: ignore
                    ):
                        states[str(mate_1.query_name)].append(
                            StatesOrientationRelationship.
                            inward_stranded_reverse
                        )
        except (OSError, ValueError) as exc:
            raise FileProblem(
                f"Failed to open SAM file: '{sam}'"
            ) from exc

        # deciding read orientation
        LOGGER.debug("Deciding read orientation...")
        reads = len(states)
        fractions = [
            self.get_frequencies(*state) for state in states.values()
        ]
        fractions_all_states = {
            key: val / reads for (key, val) in
            self.sum_dicts(*fractions).items()
        }
        orientation: ResultsOrientation = ResultsOrientation()
        if len(fractions_all_states) == 0:
            return orientation
        most_common_state = max(
            fractions_all_states,
            key=fractions_all_states.get  # type: ignore
        )
        fraction_most_common_state = max(fractions_all_states.values())
        if reads >= self.min_mapped_reads:
            if fraction_most_common_state > self.min_fraction:
                orientation.relationship = most_common_state
                assert orientation.relationship is not None
                if orientation.relationship.value[-2:] == "SF":
                    orientation.file_1 = StatesOrientation("SF")
                    orientation.file_2 = StatesOrientation("SR")
                else:
                    orientation.file_1 = StatesOrientation("SR")
                    orientation.file_2 = StatesOrientation("SF")
            else:
                orientation.relationship = (
                    StatesOrientationRelationship.inward_unstranded
                )
                orientation.file_1 = StatesOrientation.unstranded
                orientation.file_2 = StatesOrientation.unstranded

        orient_df_1 = self.create_orient_df(
            reads, fractions_all_states, orientation, paired=True, file_index=1
        )
        orient_df_2 = self.create_orient_df(
            reads, fractions_all_states, orientation, paired=True, file_index=2
        )

        LOGGER.debug(
            f"Required number of mapped reads: {self.min_mapped_reads}"
        )
        LOGGER.debug(f"Number of mapped reads: {orient_df_1.iloc[0, 0]}")
        LOGGER.debug(f"Fraction of ISF: {orient_df_1.iloc[0, 1]}")
        LOGGER.debug(f"Fraction of ISR: {orient_df_1.iloc[0, 2]}")
        LOGGER.debug(f"Orientation file 1: {orient_df_1.iloc[0, 3]}")
        LOGGER.debug(f"Orientation file 2: {orient_df_2.iloc[0, 3]}")
        LOGGER.debug(
            f"Orientation relationship: {orient_df_1.iloc[0, 4]}"
        )

        self.write_orientation_to_json(
            orient_df_1, getattr(self.paths[0], 'name')
        )
        self.write_orientation_to_json(
            orient_df_2, getattr(self.paths[1], 'name')
        )

        return orientation

    @staticmethod
    def get_frequencies(*items: Any) -> Dict[Any, float]:
        """Get frequencies of arguments as fractions of the number of all
        arguments.

        Args:
            *items: Items to get frequencies for.

        Returns:
            Dictionary of arguments and their frequencies.
        """
        counts: DefaultDict[Any, int] = defaultdict(int)
        length: int = len(items)
        for item in items:
            counts[item] += 1
        fractions: Dict[Any, float] = {
            k: v/length for (k, v) in counts.items()
        }
        return fractions

    @staticmethod
    def sum_dicts(*dicts: Dict[Any, float]) -> Dict[Any, float]:
        """Sum of dictionaries with numeric values.

        Args:
            *dicts: Dictionaries to sum up.

        Returns:
            Dictionary with union of keys of input dictionaries and all values
            added up.
        """
        result: DefaultDict[Any, float] = defaultdict(int)
        for dct in dicts:
            for key, num in dct.items():
                result[key] += num
        return dict(result)

    @staticmethod
    def create_orient_df(
            reads,
            fractions_all_states,
            orientation,
            paired: bool,
            file_index=None
    ):
        """Prepare DataFrame for orientation details.

        Constructs a DataFrame with information about read orientation for
        single or paired-end sequencing data.

        Args:
            reads: Number of mapped reads.
            fractions_all_states: Dictionary containing the fraction
                of each orientation state.
            orientation: Orientation states.
            paired: Indicates if the sequencing data is paired-end.
            file_index: Specifies the index of the file for paired-end data
                (1 or 2). Ignored for single-end data.

        Returns:
            pd.DataFrame: A DataFrame containing orientation details.
        """
        if paired:
            data = {
                'Number of mapped reads': reads,
                'Fraction ISF': fractions_all_states.get(
                    StatesOrientationRelationship.inward_stranded_forward
                ),
                'Fraction ISR': fractions_all_states.get(
                    StatesOrientationRelationship.inward_stranded_reverse
                ),
                'Orientation': getattr(
                    orientation.file_1
                    if file_index == 1 else orientation.file_2,
                    'value',
                    None
                ),
                'Relationship': getattr(
                    orientation.relationship, 'value', None
                )
            }
        else:
            data = {
                'Number of mapped reads': reads,
                'Fraction SF': fractions_all_states.get(
                    StatesOrientation.stranded_forward
                ),
                'Fraction SR': fractions_all_states.get(
                    StatesOrientation.stranded_reverse
                ),
                'Orientation': orientation.value
            }
        return pd.DataFrame([data])

    def write_orientation_to_json(self, orient_df, filename):
        """Write orientation dataframe to a JSON file.

        Serializes the provided orientation dataframe to a JSON file
            with indentation.

        Args:
            orient_df: The dataframe containing orientation details.
            filename: Name of the file to save the JSON data.

        Returns:
            None
        """
        file_path = Path(self.out_dir) / f"read_orientation_{filename}.json"
        LOGGER.debug(f"Writing results to file: {file_path}")
        orient_df.to_json(file_path, orient='split', index=False, indent=True)
