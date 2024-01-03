"""Infer mate information from sample data."""

import logging
from pathlib import Path
from typing import (List, Optional)
import re
import pysam

from Bio.SeqIO.QualityIO import FastqGeneralIterator  # type: ignore

from htsinfer.exceptions import (
    FileProblem,
    InconsistentFastqIdentifiers,
    MetadataWarning,
    UnknownFastqIdentifier,
)
from htsinfer.models import (
    ResultsType,
    StatesType,
    StatesTypeRelationship,
    SeqIdFormats,
    Config,
)
from htsinfer.mapping import Mapping

LOGGER = logging.getLogger(__name__)


class GetLibType:
    """Determine type (single/paired) information for a single or a pair of
    FASTQ sequencing libraries.

    Args:
        config: Container class for all arguments used in inference
                and results produced by the class.

    Attributes:
        path_1: Path to single-end library or first mate file.
        path_2: Path to second mate file.
        results: Results container for storing library type information for
            the provided files, as well as the mate relationship between the
            two files, if applicable.

    Examples:
        >>> GetLibType(
        ...     path_1="tests/files/first_mate.fastq"
        ... ).evaluate()
        ResultsType(file_1=<OutcomesType.single: 'single'>, file_2=<OutcomesTyp
e.not_available: 'not_available'>, relationship=<OutcomesTypeRelationship.not_a
vailable: 'not_available'>)
        >>> GetLibType(
        ...     path_1="tests/files/first_mate.fastq",
        ...     path_2="../tests/test_files/second_mate.fastq",
        ... ).evaluate()
        ResultsType(file_1=<OutcomesType.first_mate: 'first_mate'>, file_2=<Out
comesType.second_mate: 'second_mate'>, relationship=<OutcomesTypeRelationship.s
plit_mates: 'split_mates'>)
        ('first_mate', 'second_mate', 'split_mates')
    """
    def __init__(
        self,
        config: Config,
        mapping: Mapping,
    ):
        """Class constructor."""
        self.path_1: Path = config.args.path_1_processed
        self.path_2: Optional[Path] = config.args.path_2_processed
        self.library_source = config.results.library_source
        self.results: ResultsType = ResultsType()
        self.tmp_dir = config.args.tmp_dir
        self.mapping = mapping
        self.max_distance = config.args.lib_type_max_distance
        self.cutoff = config.args.lib_type_mates_cutoff

    def evaluate(self) -> None:
        """Decide type information and mate relationship."""

        # process file 1
        LOGGER.debug(f"Processing file: '{self.path_1}'")
        mates_file_1 = GetFastqType(path=self.path_1)
        mates_file_1.evaluate()
        self.results.file_1 = mates_file_1.result
        LOGGER.debug(f"Library type: {self.results.file_1}")

        # process file 2
        if self.path_2 is not None:
            LOGGER.debug(f"Processing putative mate file: '{self.path_2}'")
            mates_file_2 = GetFastqType(path=self.path_2)
            mates_file_2.evaluate()
            self.results.file_2 = mates_file_2.result
            LOGGER.debug(f"Library type: {self.results.file_2}")

            # check whether libraries are from a pair
            LOGGER.debug("Checking mate relationship between files...")
            self._evaluate_mate_relationship(
                ids_1=mates_file_1.seq_ids,
                ids_2=mates_file_2.seq_ids,
            )
            LOGGER.debug(f"Mate relationship: {self.results.relationship}")

    def _evaluate_mate_relationship(
        self,
        ids_1: List[str],
        ids_2: List[str],
    ) -> None:
        """Decide mate relationship.

        Args:
            ids_1: List of sequence identifier prefixes of the putative first
            mate file, i.e., the fragments up until the mate information,
            if available, as defined by a named capture group ``prefix`` in a
            regular expression to extract mate information.
            ids_2: As `ids_1` but for the putative second mate file.
        """
        self.results.relationship = StatesTypeRelationship.not_mates
        if ids_1 == ids_2:
            if (
                self.results.file_1 == StatesType.first_mate and
                self.results.file_2 == StatesType.second_mate
            ) or (
                self.results.file_1 == StatesType.second_mate and
                self.results.file_2 == StatesType.first_mate
            ):
                self.results.relationship = (
                    StatesTypeRelationship.split_mates
                )
                self.mapping.library_type.relationship = (
                    StatesTypeRelationship.split_mates
                )
        else:
            self.mapping.library_type.relationship \
                = StatesTypeRelationship.not_available
            self.mapping.library_source = self.library_source
            self.mapping.paths = self.path_1, self.path_2
            self.mapping.evaluate()
            self._align_mates()

    def _align_mates(self):
        """Decide mate relationship by alignment."""

        alignment_1 = self.mapping.star_dirs[0] / 'Aligned.out.sam'
        alignment_2 = self.mapping.star_dirs[1] / 'Aligned.out.sam'

        samfile1 = pysam.AlignmentFile(str(alignment_1), 'r')
        samfile2 = pysam.AlignmentFile(str(alignment_2), 'r')

        seq_id1 = None
        seq_id2 = None
        previous_seq_id1 = None
        previous_seq_id2 = None

        reads1 = []  # List to store alignments for one read from file1
        mate1 = []  # List to store alignments for each read
        reads2 = []  # List to store alignments for one read from file2

        concordant = 0

        for read1 in samfile1:
            if not read1.flag & (1 << 2) and isinstance(read1.query_name, str):
                seq_id1 = read1.query_name
                if (
                    seq_id1 != previous_seq_id1 and
                    previous_seq_id1 is not None
                ):
                    mate1.append(reads1.copy())
                    reads1.clear()
                if read1.reference_end:
                    reads1.append(read1)
            previous_seq_id1 = seq_id1
        mate1.append(reads1.copy())

        read_counter = 0
        for read2 in samfile2:
            if not read2.flag & (1 << 2) and isinstance(read2.query_name, str):
                seq_id2 = read2.query_name
                if seq_id2 != previous_seq_id2 \
                        and previous_seq_id2 is not None:
                    if self._compare_alignments(mate1[read_counter], reads2):
                        concordant += 1
                    reads2.clear()
                    read_counter += 1
                if read2.reference_end:
                    reads2.append(read2)
            previous_seq_id2 = seq_id2

        if self._compare_alignments(mate1[read_counter], reads2):
            concordant += 1

        self._update_relationship(concordant, read_counter)

        samfile1.close()
        samfile2.close()

    def _update_relationship(self, concordant, read_counter):
        """Helper function to update relationship based on alignment."""
        try:
            if (concordant / read_counter) >= self.cutoff:
                self.results.relationship = (
                    StatesTypeRelationship.split_mates
                )
                self.mapping.library_type.relationship \
                    = StatesTypeRelationship.split_mates
                self.mapping.mapped = False
                self.mapping.star_dirs = []
            else:
                self.results.relationship = (
                    StatesTypeRelationship.not_mates
                )
        except ZeroDivisionError:
            self.results.relationship = (
                StatesTypeRelationship.not_available
            )

    class AlignedSegment:
        """Placeholder class for mypy "Missing attribute"
        error in _compare_alignments(), the actual object used
        is pysam.AlignedSegment class.
        """
        def __init__(self):
            """Class constructor."""
            self.reference_name: str = ""
            self.reference_start: int = 0
            self.reference_end: int = 0

    def _compare_alignments(
            self,
            read1: List[AlignedSegment],
            read2: List[AlignedSegment],  # List[object],
    ) -> bool:
        """Compare two lists of alignments and decide
        if they come from a single fragment.

        Args:
        read1: List of pysam.AlignedSegment records for first mate
        read2: List of pysam.AlignedSegment records for second mate
        """
        for one in read1:
            for two in read2:
                if one.reference_name == two.reference_name:
                    second_start = one.reference_start \
                        if one.reference_start > two.reference_start \
                        else two.reference_start
                    first_end = one.reference_end \
                        if one.reference_end < two.reference_end \
                        else two.reference_end
                    if second_start - first_end < self.max_distance:
                        return True
        return False


class GetFastqType():
    """Determine type (single/paired) information for an individual FASTQ
    sequencing library.

    Args:
        path: File path to read library.

    Attributes:
        path: File path to read library.
        seq_ids: List of sequence identifier prefixes of the provided read
            library, i.e., the fragments up until the mate information,
            if available, as defined by a named capture group ``prefix`` in a
            regular expression to extract mate information.
        seq_id_format: The sequence identifier format of the read library, as
            identified by inspecting the first read and matching one of the
            available regular expressions for the different identifier formats.
        result: The current best guess for the type of the provided library.

    Examples:
        >>> lib_type = GetFastqType(
        ...     path="tests/files/first_mate.fastq"
        ... ).evaluate()
        <OutcomesType.first_mate: 'first_mate'>
    """
    def __init__(
        self,
        path: Path,
    ):
        """Class constructor."""
        self.path: Path = path
        self.seq_ids: List[str] = []
        self.seq_id_format: Optional[SeqIdFormats] = None
        self.result: StatesType = StatesType.not_available

    def evaluate(self) -> None:
        """Decide library type.

        Raises:
            NoMetadataDetermined: Type information could not be determined.
        """
        records: int = 0

        try:
            with open(self.path, encoding="utf-8") as _f:  # type: ignore

                # Get sequence identifier format from first record
                LOGGER.debug(
                    "Determining identifier and library type from first "
                    "record..."
                )
                try:
                    seq_iter = FastqGeneralIterator(source=_f)
                    seq_id = next(seq_iter)[0]
                    records += 1
                    for seq_id_format in SeqIdFormats:
                        try:
                            self._get_read_type(
                                seq_id=seq_id,
                                regex=seq_id_format.value,
                            )
                        except UnknownFastqIdentifier:
                            continue
                        self.seq_id_format = seq_id_format
                        break

                except StopIteration as exc:
                    self.result = StatesType.not_available
                    raise FileProblem(f"File is empty: {self.path}") from exc

                if self.seq_id_format is None:
                    self.result = StatesType.not_available
                    raise MetadataWarning(
                        "Could not determine sequence identifier format."
                    )
                LOGGER.debug(
                    f"Sequence identifier format: {self.seq_id_format.name}"
                )

                # Ensure that remaining records are compatible with sequence
                # identifier format and library type determined from first
                # record
                LOGGER.debug(
                    "Checking consistency of remaining reads with initially "
                    "determined identifier format and library type..."
                )
                for record in seq_iter:
                    records += 1
                    try:
                        self._get_read_type(
                            seq_id=record[0],
                            regex=self.seq_id_format.value,
                        )
                    except (
                        InconsistentFastqIdentifiers,
                        UnknownFastqIdentifier,
                    ) as exc:
                        self.result = StatesType.not_available
                        raise MetadataWarning(
                            f"{type(exc).__name__}: {str(exc)}"
                        ) from exc

        except (OSError, ValueError) as exc:
            self.result = StatesType.not_available
            raise FileProblem(f"{type(exc).__name__}: {str(exc)}") from exc

        LOGGER.debug(f"Total records processed: {records}")

    def _get_read_type(
        self,
        seq_id: str,
        regex: re.Pattern,
    ) -> None:
        """Get/update library type information from sequence identifier.

        Args:
            seq_id: Sequence identifier.
            regex: A regular expression to extract library type information
                from a read. The expression needs to contain a named group
                ``mate`` that MAY be present in the sequence identifier and
                that, if present, MUST match a single character that can take
                values ``1`` (for first-mate files) and ``2`` (for second-mate
                files), as well as a named group ``prefix`` that MUST be
                present in the sequence identifier and that contains all of the
                sequence identifier that is identical between a mate pair up
                until the mate information itself.

        Raises:
            InconsistentFastqIdentifiers: A sequence identifier was encountered
                that suggests a different library type than previous
                identifiers.
            UnkwownFastqIdentifier: A sequence identifier of unknown format was
                encountered.
        """
        # Note: Conditionals have been optimized to minimize required checks
        # for the most likely scenarios, not to simplify code
        match = re.search(regex, seq_id)
        if match:
            self.seq_ids.append(match.group('prefix'))
            # Read appears to be derived from a single-end library
            if match.group('mate') is None:
                self._get_read_type_single(seq_id=seq_id)
            # Read appears to be derived from a paired-end library
            else:
                # First mate
                if int(match.group('mate')) == 1:
                    self._get_read_type_paired_mate_1(seq_id=seq_id)
                # Second mate
                else:
                    self._get_read_type_paired_mate_2(seq_id=seq_id)
        else:
            raise UnknownFastqIdentifier(
                f"Cannot determine identifier format: {seq_id}"
            )

    def _get_read_type_single(
        self,
        seq_id: str,
    ) -> None:
        """Helper function to process single-end libraries in
        ``GetFastqType._get_read_type()``.

        Args:
            seq_id: Sequence identifier.

        Raises:
            InconsistentFastqIdentifiers: A sequence identifier was encountered
                that suggests a different library type than previous
                identifiers.
        """
        if self.result == StatesType.single:
            pass
        elif self.result == StatesType.not_available:
            self.result = StatesType.single
        else:
            raise InconsistentFastqIdentifiers(
                "Previous sequence identifiers suggest that reads are part of "
                "a paired-end sequencing library, but current identifier "
                f"suggests a single-end library: {seq_id}"
            )

    def _get_read_type_paired_mate_1(
        self,
        seq_id: str,
    ) -> None:
        """Helper function to process first mate files of paired-end libraries
        in ``GetFastqType._get_read_type()``.

        Args:
            seq_id: Sequence identifier.

        Raises:
            InconsistentFastqIdentifiers: A sequence identifier was encountered
                that suggests a different library type than previous
                identifiers.
        """
        if self.result in (
            StatesType.first_mate,
            StatesType.mixed_mates
        ):
            pass
        elif self.result == StatesType.second_mate:
            self.result = StatesType.mixed_mates
        elif self.result == StatesType.not_available:
            self.result = StatesType.first_mate
        else:
            raise InconsistentFastqIdentifiers(
                "Previous sequence identifiers suggest that reads are part of "
                "a single-end sequencing library, but current identifier "
                f"suggests a paired-end library: {seq_id}"
            )

    def _get_read_type_paired_mate_2(
        self,
        seq_id: str,
    ) -> None:
        """Helper function to process second mate files of paired-end libraries
        in ``GetFastqType._get_read_type()``.

        Args:
            seq_id: Sequence identifier.

        Raises:
            InconsistentFastqIdentifiers: A sequence identifier was encountered
                that suggests a different library type than previous
                identifiers.
        """
        if self.result in (
            StatesType.second_mate,
            StatesType.mixed_mates
        ):
            pass
        elif self.result == StatesType.first_mate:
            self.result = StatesType.mixed_mates
        elif self.result == StatesType.not_available:
            self.result = StatesType.second_mate
        else:
            raise InconsistentFastqIdentifiers(
                "Previous sequence identifiers suggest that reads are part of "
                "a single-end sequencing library, but current identifier "
                f"suggests a paired-end library: {seq_id}"
            )
