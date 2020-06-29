"""Unit tests for module 'htsinfer.py'"""

import importlib.util
import logging
import os

import pytest

from htsinfer.htsinfer import (
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
    "htsinfer",
    "htsinfer.py",
)
OPTION_HELP = "help"
OPTION_SWITCH = "debug"
OPTION_SWITCH_DEFAULT = False
OPTION_INVALID = "-+-+-"
USER_INPUT = "user_input"
MOCK_FILE_PATH = "/path/to/file"
TEST_FILE = os.path.join(
    os.path.dirname(__file__),
    "sample_files",
    "first_mate.fastq",
)


# main()
def test_main_as_script():
    fl = os.path.join(os.path.dirname(__file__), MAIN_FILE)
    spec = importlib.util.spec_from_file_location('__main__', fl)
    module = importlib.util.module_from_spec(spec)
    with pytest.raises(SystemExit):
        spec.loader.exec_module(module)


def test_main_with_args():
    args = parse_args([
        "--file-1", TEST_FILE,
        "--file-2", TEST_FILE,
        "--debug",
    ])
    ret = main(args=args)
    assert ret is None


# parse_args()
def test_help_option():
    with pytest.raises(SystemExit):
        assert parse_args(["--" + OPTION_HELP])


def test_no_args():
    with pytest.raises(SystemExit):
        parse_args([])


def test_switch():
    ret = parse_args([
        "--file-1", MOCK_FILE_PATH,
        "--" + OPTION_SWITCH,
    ])
    assert vars(ret)[OPTION_SWITCH] is not OPTION_SWITCH_DEFAULT


def test_invalid_option():
    with pytest.raises(SystemExit):
        assert parse_args([OPTION_INVALID])


# setup_logging()
def test_log_level_default():
    logger = logging.getLogger("my_logger")
    setup_logging()


def test_log_level_verbose():
    logger = logging.getLogger("my_logger")
    setup_logging(
        verbose=True,
    )


def test_log_level_debug():
    logger = logging.getLogger("my_logger")
    setup_logging(
        debug=True,
    )

