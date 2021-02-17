#!/usr/bin/env python
"""HTSinfer infers metadata from High Throughput Sequencing (HTS) data"""

__version__ = "0.1.0"
__copyright__ = "Copyright 2020 Zavolan lab, Biozentrum, University of Basel"
__license__ = "Apache license 2.0"
__author__ = "Rohan Kandhari"
__maintainer__ = "Rohan Kandhari"
__email__ = "rohan.kandhari.bme16@iitbhu.ac.in"

import argparse
import logging
import sys
from typing import (Optional, Sequence)

from htsinfer import (
    infer_single_paired,
    infer_organism,
    infer_adapter
)

logger = logging.getLogger(__name__)


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
        '-mm', '--min-match',
        metavar="FLOAT",
        type=float,
        default=10,
        help=(
            "minimum match percentage that organism needs to have"
        )
    )
    parser.add_argument(
        '-fr', '--factor',
        metavar="FLOAT",
        type=float,
        default=2,
        help=(
            "factor by which first organism is greater than the second"
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
    args = parse_args()
    setup_logging(
        verbose=args.verbose,
        debug=args.debug,
    )
    logger.info("Started script...")
    logger.debug(f"CLI options: {args}")
    results = {}
    results['single_paired'] = infer_single_paired.infer(
        file_1=args.file_1,
        file_2=args.file_2,
    )
    '''
    results['organism'] = infer_organism.kallisto(
        file_1=args.file_1,
        file_2=args.file_2,
        min_match=args.min_match,
        factor=args.factor
    )
    '''
    results['adapters'] = infer_adapter.infer(
        file_1=args.file_1,
        file_2=args.file_2
    )
    logger.info(f"Results: {results}")
    logger.info("Done.")


if __name__ == "__main__":
    main()
