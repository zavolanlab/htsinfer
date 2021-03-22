#!/usr/bin/env python
"""Command-line interface client."""

import argparse
from enum import Enum
import logging
from pathlib import Path
import signal
import sys

from htsinfer import (HtsInfer, __version__)

LOGGER = logging.getLogger(__name__)


class LogLevels(Enum):
    """Log level enumerator."""
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARN = logging.WARNING
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments.

    Returns:
        Parsed CLI arguments.
    """
    # set metadata
    usage = (
        """htsinfer [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}] [-h]
                [--version]
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
