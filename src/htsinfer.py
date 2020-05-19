#!/usr/bin/env python
"""HTSinfer infers metadata from High Throughput Sequencing (HTS) data"""

__version__ = "0.1.0"
__copyright__ = "Copyright 2020 Zavolan lab, Biozentrum, University of Basel"
__license__ = "Apache license 2.0"
__author__ = "Rohan Kandhari"
__maintainer__ = "Rohan Kandhari"
__email__ = "rohan.kandhari.bme16@iitbhu.ac.in"

# TODO AUTHOR: add here built-in modules
import argparse
import logging
import sys
from typing import (Optional, Sequence)

# TODO AUTHOR: add here third party modules

# TODO AUTHOR: add here own modules

LOGGER = logging.getLogger(__name__)


def parse_args(
    args: Optional[Sequence[str]] = None
) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        # TODO AUTHOR: add here detailed tool description; leave a blank line
        # in between synopsis and extended description
        description=sys.modules[__name__].__doc__,
    )

    # TODO AUTHOR: add here optional and positional arguments as per argparse
    # docs; for many optional arguments, consider adding argument groups for
    # clarity
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
        help="also print debugging messages to STDERR",
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
        help="show version information and exit",
    )

    return parser.parse_args(args)


def setup_logging(
    logger: logging.Logger = logging.getLogger(__name__),
    verbose: bool = False,
    debug: bool = False,
) -> logging.Logger:
    """Configure logging."""
    if debug:
        logger.setLevel(logging.DEBUG)
    elif verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(
        "[%(asctime)-15s: %(levelname)-8s @ %(funcName)s] %(message)s"
    ))
    logger.addHandler(handler)

    return logger


def main():
    """Main function."""
    args = parse_args()
    logger = setup_logging(
        logger=LOGGER,
        verbose=args.verbose,
        debug=args.debug,
    )
    # TODO AUTHOR: Put main code here. Options and positional arguments are in
    # `args`, logging can be used with `logger`; see (and delete) example log
    # message below
    logger.info("Started script.")
    logger.info(f"CLI options: {args}")


if __name__ == "__main__":
    main()
