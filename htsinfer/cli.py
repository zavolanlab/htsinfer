#!/usr/bin/env python
"""Command-line interface client."""

import argparse
import logging
from pathlib import Path
import signal
import sys
import tempfile

from htsinfer import (
    HtsInfer,
    __version__,
)
from htsinfer.models import (
    CleanupRegimes,
    LogLevels,
    Args,
    Results,
    Config,
)

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments.

    Returns:
        Parsed CLI arguments.
    """
    # set metadata
    description = (
        f"{sys.modules[__name__].__doc__}\n\n"
        ""
    )
    epilog = (
        f"%(prog)s v{__version__}, (c) 2021 by Zavolab "
        "(zavolab-biozentrum@unibas.ch)"
    )

    # custom actions
    class PathsAction(argparse.Action):
        """Sanitize ``paths`` parsing in positional args."""
        def __call__(
            self,
            parser,
            namespace,
            values,
            option_string=None,
        ) -> None:
            if len(values) > 2:
                parser.print_usage(file=sys.stderr)
                sys.stderr.write(
                    "htsinfer: error: only one or two of the following "
                    "arguments are allowed: PATH\n"
                )
                parser.exit(2)
            if len(values) == 1:
                values.append(None)
            setattr(namespace, self.dest, values)

    # instantiate parser
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # add arguments
    parser.add_argument(
        'paths',
        nargs="+",
        type=Path,
        action=PathsAction,
        metavar="PATH",
        help=(
            "either one or two paths to FASTQ files representing the "
            "sequencing library to be evaluated, for single- or paired-ended "
            "libraries, respectively"
        ),
    )
    parser.add_argument(
        "--output-directory",
        default=Path.cwd() / 'results_htsinfer',
        type=lambda p: Path(p).absolute(),
        metavar="PATH",
        help="path to directory where output is written to",
    )
    parser.add_argument(
        "--temporary-directory",
        default=Path(tempfile.gettempdir()) / 'tmp_htsinfer',
        type=lambda p: Path(p).absolute(),
        metavar="PATH",
        help="path to directory where temporary output is written to",
    )
    parser.add_argument(
        "--cleanup-regime",
        choices=[e.name for e in CleanupRegimes],
        default=CleanupRegimes.DEFAULT.name,
        type=str,
        help=(
            "determine which data to keep after each run; in default mode, "
            "both temporary data and results are kept when '--verbosity' is "
            "set to 'DEBUG', no data is kept when all metadata could be "
            "successfully determined, and only results are kept otherwise"
        )
    )
    parser.add_argument(
        "--records",
        default=0,
        type=int,
        metavar="INT",
        help=(
            "number of records to process; if set to ``0`` or if the "
            "specified value equals or exceeds the number of available "
            "records, all records will be processed"
        )
    )
    parser.add_argument(
        "--threads",
        metavar="INT",
        type=int,
        default=1,
        help=(
            "number of threads to run STAR with"
        )
    )
    parser.add_argument(
        "--transcripts",
        metavar="FASTA",
        type=str,
        default=(
            Path(__file__).parents[1].absolute() / "data/transcripts.fasta.gz"
        ),
        help=(
            "FASTA file containing transcripts to be used for mapping files "
            "`--file-1` and `--file-2` for inferring library source and read "
            "orientation. Requires that sequence identifier lines are "
            "separated by the pipe (`|`) character and that the 4th and 5th "
            "columns contain a short source name and taxon identifier, "
            "respectively. Example sequence identifier: "
            "`rpl-13|ACYPI006272|ACYPI006272-RA|apisum|7029`"
        )
    )
    parser.add_argument(
        "--read-layout-adapters",
        metavar="PATH",
        type=str,
        default=(
            Path(__file__).parent.parent.absolute() /
            "data/adapter_fragments.txt"
        ),
        help=(
            "path to text file containing 3' adapter sequences to scan for "
            "(one sequence per line)"
        )
    )
    parser.add_argument(
        "--read-layout-min-match-percentage",
        metavar="FLOAT",
        type=float,
        default=2,
        help=(
            "minimum percentage of reads that contain a given adapter "
            "sequence in order for it to be considered as the library's "
            "3'-end adapter"
        )
    )
    parser.add_argument(
        "--read-layout-min-frequency-ratio",
        metavar="FLOAT",
        type=float,
        default=2,
        help=(
            "minimum frequency ratio between the first and second most "
            "frequent adapter in order for the former to be considered as the "
            "library's 3'-end adapter"
        )
    )
    parser.add_argument(
        "--library-source-min-match-percentage",
        metavar="FLOAT",
        type=float,
        default=5,
        help=(
            "minimum percentage of reads that are consistent with a given "
            "source in order for it to be considered the library's source"
        )
    )
    parser.add_argument(
        "--library-source-min-frequency-ratio",
        metavar="FLOAT",
        type=float,
        default=2,
        help=(
            "minimum frequency ratio between the first and second "
            "most frequent source in order for the former to be considered "
            "the library's source"
        )
    )
    parser.add_argument(
        "--read-orientation-min-mapped-reads",
        metavar="INT",
        type=int,
        default=20,
        help=(
            "minimum number of mapped reads for deeming the read orientation "
            "result reliable"
        )
    )
    parser.add_argument(
        "--read-orientation-min-fraction",
        metavar="FLOAT",
        type=float,
        default=0.75,
        help=(
            "minimum fraction of mapped reads required to be consistent with "
            "a given read orientation state in order for that orientation to "
            "be reported. Must be above 0.5."
        )
    )
    parser.add_argument(
        "--verbosity",
        choices=[e.name for e in LogLevels],
        default=LogLevels.INFO.name,
        type=str,
        help="logging verbosity level",
    )
    parser.add_argument(
        "-h", "--help",
        action="help",
        help="show this help message and exit",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=epilog,
        help="show version information and exit",
    )

    # fix faulty usage string for nargs={1,2}
    parser.usage = parser.format_usage()
    parser.usage = parser.usage.replace("PATH [PATH ...]", "PATH [PATH]")

    # return parsed arguments
    return parser.parse_args()


def setup_logging(verbosity: str = 'INFO') -> None:
    """Configure logging.

    Args:
        verbosity: Level of logging verbosity.
    """
    level = LogLevels[verbosity].value
    logging.basicConfig(
        level=level,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt='%Y-%m-%d %H:%M:%S',
    )


def main() -> None:
    """Entry point for CLI executable."""
    try:
        # handle CLI args
        args = parse_args()

        arguments = Args(
            path_1=args.paths[0],
            path_2=args.paths[1],
            out_dir=args.output_directory,
            tmp_dir=args.temporary_directory,
            cleanup_regime=CleanupRegimes[args.cleanup_regime],
            records=args.records,
            threads=args.threads,
            transcripts_file=args.transcripts,
            read_layout_adapter_file=args.read_layout_adapters,
            read_layout_min_match_pct=args.read_layout_min_match_percentage,
            read_layout_min_freq_ratio=args.read_layout_min_frequency_ratio,
            lib_source_min_match_pct=args.library_source_min_match_percentage,
            lib_source_min_freq_ratio=args.library_source_min_frequency_ratio,
            read_orientation_min_mapped_reads=(
                args.read_orientation_min_mapped_reads
            ),
            read_orientation_min_fraction=(
                args.read_orientation_min_fraction
            ),
        )

        results = Results()

        config = Config(
            args=arguments,
            results=results,
        )

        # set up logging
        setup_logging(verbosity=args.verbosity)
        LOGGER.info("Started HTSinfer...")
        LOGGER.debug(f"CLI arguments: {args}")

        # determine library metadata
        hts_infer = HtsInfer(
            config=config,
        )
        hts_infer.evaluate()
        hts_infer.print()

    except KeyboardInterrupt:
        LOGGER.error('Execution interrupted.')
        sys.exit(128 + signal.SIGINT)

    # conclude execution
    LOGGER.info("Done")
    sys.exit(hts_infer.state.value)


if __name__ == '__main__':
    main()
