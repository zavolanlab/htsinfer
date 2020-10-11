"""Unit tests for module 'htsinfer.py'"""

import importlib.util
import logging
import os
import sys

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
    "test_files",
    "first_mate.fastq",
)


def test_main_as_script():
    """Run as script."""
    fl = os.path.join(os.path.dirname(__file__), MAIN_FILE)
    spec = importlib.util.spec_from_file_location('__main__', fl)
    module = importlib.util.module_from_spec(spec)
    with pytest.raises(SystemExit):
        spec.loader.exec_module(module)


class TestMain:
    """Test `main()` function."""

    def test_main(self):
        """Call without args."""
        with pytest.raises(SystemExit):
            main()

    def test_main_with_args(self, monkeypatch):
        """Call with args."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                '--file-1', TEST_FILE,
                '--file-2', TEST_FILE,
            ]
        )
        assert main() is None


class TestParseArgs:
    """Test `parse_args()` function."""

    def test_help_option(self):
        """Test help option."""
        with pytest.raises(SystemExit):
            assert parse_args(["--" + OPTION_HELP])

    def test_no_args(self):
        """Call without args."""
        with pytest.raises(SystemExit):
            parse_args([])

    def test_switch(self):
        """Test switch parameter."""
        ret = parse_args([
            "--file-1", MOCK_FILE_PATH,
            "--" + OPTION_SWITCH,
        ])
        assert vars(ret)[OPTION_SWITCH] is not OPTION_SWITCH_DEFAULT

    def test_invalid_option(self):
        """Test invalid argument."""
        with pytest.raises(SystemExit):
            assert parse_args([OPTION_INVALID])


class TestSetupLogging:
    """Test `setup_logging()` function."""

    def test_log_level_default(self):
        "Call without args."
        setup_logging()

    def test_log_level_verbose(self):
        "Set log level `logging.INFO`."
        setup_logging(
            verbose=True,
        )

    def test_log_level_debug(self):
        "Set log level `logging.DEBUG`."
        setup_logging(
            debug=True,
        )
