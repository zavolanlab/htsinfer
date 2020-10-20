#!/usr/bin/env python
"""Infer experiment metadata from High Throughput Sequencing (HTS) data."""

import argparse
import logging
import sys
from typing import (Any, Dict, Optional, Sequence)

from htsinfer import (
    infer_single_paired,
    infer_read_layout,
    __version__,
)

LOGGER = logging.getLogger(__name__)


def parse_args(
    args: Optional[Sequence[str]] = None
) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
    )

    parser.add_argument(
        '-f1', '--file-1',
        metavar="FILE",
        type=str,
        required=True,
        help="file path to read/first mate library",
    )
    parser.add_argument(
        '-f2', '--file-2',
        metavar="FILE",
        type=str,
        default=None,
        help="file path to second mate library",
    )
    parser.add_argument(
        '-n', '--max-records',
        metavar="INT",
        type=int,
        default=10000,
        help=(
            "maximum number of records to process, starting with first "
            "record; set to 0 to process entire file(s)"
        )
    )
    parser.add_argument(
        '-o', '--organism',
        metavar="STR",
        type=str,
        default='hsapiens',
        help=(
            "source organism of the sequencing library, either as a short "
            "name (e.g., 'hsapiens') or taxon identifier (e.g., '9606'); if "
            "provided, will not not be inferred by the application"
        )
    )
    parser.add_argument(
        '--verbose', "-v",
        action='store_true',
        default=False,
        help="print logging messages to STDERR",
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        default=False,
        help="print debugging messages to STDERR",
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
        help="show version information and exit",
    )

    return parser.parse_args(args)


def setup_logging(
    verbose: bool = False,
    debug: bool = False,
) -> None:
    """Configure logging."""
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(message)s",
        datefmt='%m-%d %H:%M:%S',
    )


def main() -> None:
    """Main function.

    Args:
        args: Command-line arguments and their values.
    """

    # Parse CLI arguments
    args = parse_args()

    # Set up logging
    setup_logging(
        verbose=args.verbose,
        debug=args.debug,
    )

    # Initialize results container
    LOGGER.info("Started script...")
    LOGGER.debug(f"CLI options: {args}")
    results: Dict[str, Any] = {}

    # Infer library type
    results['single_paired'] = infer_single_paired.infer(
        file_1=args.file_1,
        file_2=args.file_2,
    )

    # Infer read layout
    results['read_layout'] = infer_read_layout.infer(
        file_1=args.file_1,
        file_2=args.file_2,
        organism=args.organism,
    )

    # Log results & end script
    LOGGER.info(f"Results: {results}")
    LOGGER.info("Done.")


if __name__ == "__main__":
    main()
