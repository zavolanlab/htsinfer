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
)

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments.

    Returns:
        Parsed CLI arguments.
    """
    # set metadata
    usage = (
        """htsinfer [--output-directory PATH] [--temporary-directory PATH]
                [--cleanup-regime {DEFAULT,KEEP_ALL,KEEP_NONE,KEEP_RESULTS}]
                [--records INT ] [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}]
                [-h] [--version]
                FASTQ_PATH [FASTQ_PATH]
        """
    )
    description = (
        f"{sys.modules[__name__].__doc__}\n\n"
        ""
    )
    epilog = (
        f"%(prog)s v{__version__}, (c) 2021 by Zavolab "
        "(zavolab-biozentrum@unibas.ch)"
    )

    # instantiate parser
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        usage=usage,
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # add arguments
    parser.add_argument(
        'paths',
        nargs="+",
        type=Path,
        metavar="FASTQ_PATH",
        help=(
            "either one or two file paths to the read library to be "
            "evaluated."
        ),
    )
    parser.add_argument(
        "--output-directory",
        default=Path.cwd(),
        type=lambda p: Path(p).absolute(),
        metavar="PATH",
        help="path to directory where output is written to",
    )
    parser.add_argument(
        "--temporary-directory",
        default=Path(tempfile.gettempdir()),
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

    # return parsed arguments
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments.

    Args:
        args: CLI arguments.

    Raises:
        ValueError: A CLI argument is invalid.
    """
    # ensure that not more than two file paths are passed
    if len(args.paths) > 2:
        raise ValueError("A maximum of two file paths can be specified.")
    if len(args.paths) == 1:
        args.paths.append(None)


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
        validate_args(args=args)

        # set up logging
        setup_logging(verbosity=args.verbosity)
        LOGGER.info("Started HTSinfer...")
        LOGGER.debug(f"CLI arguments: {args}")

        # determine library metadata
        hts_infer = HtsInfer(
            path_1=args.paths[0],
            path_2=args.paths[1],
            out_dir=args.output_directory,
            tmp_dir=args.temporary_directory,
            cleanup_regime=CleanupRegimes[args.cleanup_regime],
            records=args.records,
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
