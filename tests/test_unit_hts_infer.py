"""Unit tests for module 'htsinfer.py'"""

import importlib.util
import logging
import os

import pytest

from src.htsinfer import (
    main,
    parse_args,
    setup_logging,
)

# Test parameters
LOGGER = logging.getLogger(__name__)
LOG_LEVEL_DEFAULT = logging.WARNING
LOG_LEVEL_VERBOSE = logging.INFO
LOG_LEVEL_DEBUG = logging.DEBUG
MAIN_FILE = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "src",
    "htsinfer.py",
)
OPTION_HELP = "help"
OPTION_SWITCH = "debug"
OPTION_SWITCH_DEFAULT = False
OPTION_INVALID = "-+-+-"
USER_INPUT = "user_input"


# main()
def test_main():
    main()


def test_main_as_script():
    fl = os.path.join(os.path.dirname(__file__), MAIN_FILE)
    spec = importlib.util.spec_from_file_location('__main__', fl)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


# parse_args()
def test_help_option():
    with pytest.raises(SystemExit):
        assert parse_args(["--" + OPTION_HELP])


def test_no_args():
    ret = parse_args([])
    assert vars(ret)[OPTION_SWITCH] is OPTION_SWITCH_DEFAULT


def test_switch():
    ret = parse_args(["--" + OPTION_SWITCH])
    assert vars(ret)[OPTION_SWITCH] is not OPTION_SWITCH_DEFAULT


def test_invalid_option():
    with pytest.raises(SystemExit):
        assert parse_args([OPTION_INVALID])


# setup_logging()
def test_log_level_default():
    logger = setup_logging()
    assert logger.level == LOG_LEVEL_DEFAULT


def test_log_level_verbose():
    logger = setup_logging(
        verbose=True,
    )
    assert logger.level == LOG_LEVEL_VERBOSE


def test_log_level_debug():
    logger = setup_logging(
        debug=True,
    )
    assert logger.level == LOG_LEVEL_DEBUG


def test_log_level_precedence():
    logger = setup_logging(
        verbose=True,
        debug=True,
    )
    assert logger.level == LOG_LEVEL_DEBUG


def test_logger_nondefault():
    logger = setup_logging(
        logger=LOGGER,
    )
    assert logger.level == LOG_LEVEL_DEFAULT
