#!/usr/bin/env python
"""Command-line interface client."""

import argparse
import logging
from pathlib import Path
from typing import Dict
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
        dest="out_dir",
        default=Path.cwd() / 'results_htsinfer',
        type=lambda p: Path(p).absolute(),
        metavar="PATH",
        help="path to directory where output is written to",
    )
    parser.add_argument(
        "--temporary-directory",
        dest="tmp_dir",
        default=Path(tempfile.gettempdir()) / 'tmp_htsinfer',
        type=lambda p: Path(p).absolute(),
        metavar="PATH",
        help="path to directory where temporary output is written to",
    )
    parser.add_argument(
        "--cleanup-regime",
        dest="cleanup_regime",
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
        dest="records",
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
        dest="threads",
        metavar="INT",
        type=int,
        default=1,
        help=(
            "number of threads to run STAR with"
        )
    )
    parser.add_argument(
        "--transcripts",
        dest="transcripts_file",
        metavar="FASTA",
        type=str,
        default=(
            Path(__file__).parents[1].absolute()  # pylint: disable=E1101
            / "data/transcripts.fasta.gz"
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
        dest="read_layout_adapter_file",
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
        dest="read_layout_min_match_pct",
        metavar="FLOAT",
        type=float,
        default=0.1,
        help=(
            "minimum percentage of reads that contain a given adapter "
            "sequence in order for it to be considered as the library's "
            "3'-end adapter"
        )
    )
    parser.add_argument(
        "--read-layout-min-frequency-ratio",
        dest="read_layout_min_freq_ratio",
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
        dest="lib_source_min_match_pct",
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
        dest="lib_source_min_freq_ratio",
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
        dest="read_orientation_min_mapped_reads",
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
        dest="read_orientation_min_fraction",
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
        paths = args.paths
        args_dict: Dict = vars(args)
        del args_dict["paths"]
        args_dict["cleanup_regime"] = CleanupRegimes[args.cleanup_regime]

        config = Config(
            args=Args(path_1=paths[0], path_2=paths[1], **args_dict),
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
