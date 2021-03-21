"""Unit tests for module ``cli.py``."""

import argparse
import importlib.util
import os
from pathlib import Path
import sys

import pytest

from htsinfer.cli import (
    main,
    parse_args,
    setup_logging,
    validate_args,
)

# Test parameters
PACKAGE_DIR = Path(__file__).resolve().parents[1] / "htsinfer"
TEST_FILES_DIR = Path(__file__).resolve().parent / "files"
TEST_FILE = TEST_FILES_DIR / "first_mate.fastq"


class TestParseArgs:
    """Test ``parse_args()`` function."""

    def test_help_option(self, monkeypatch):
        """Test help option."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                '--help',
            ]
        )
        with pytest.raises(SystemExit) as exc:
            assert parse_args()
        assert exc.value.code == 0

    def test_no_args(self):
        """Call without args."""
        with pytest.raises(SystemExit) as exc:
            parse_args()
        assert exc.value.code == 2

    def test_positional_args(self, monkeypatch):
        """Call with positional args."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                str(TEST_FILE), str(TEST_FILE),
            ]
        )
        ret_val = parse_args()
        assert isinstance(ret_val, argparse.Namespace)


class TestValidateArgs:
    """Test ``validate_args()`` function."""

    def test_one_path(self):
        """Call with one path."""
        args = argparse.Namespace(
            paths=[TEST_FILE]
        )
        validate_args(args)
        assert len(args.paths) == 2
        assert args.paths[1] is None

    def test_two_paths(self):
        """Call with two paths."""
        args = argparse.Namespace(
            paths=[TEST_FILE, TEST_FILE]
        )
        validate_args(args)
        assert len(args.paths) == 2
        assert args.paths[1] is not None

    def test_too_many_paths(self):
        """Call with more than two paths."""
        args = argparse.Namespace(
            paths=[TEST_FILE, TEST_FILE, TEST_FILE]
        )
        with pytest.raises(ValueError):
            validate_args(args)


class TestSetupLogging:
    """Test ``setup_logging()`` function."""

    def test_log_level_default(self):
        "Call without args."
        setup_logging()

    def test_log_level_debug(self):
        "Set log level ``logging.DEBUG``."
        setup_logging(
            verbosity="DEBUG",
        )

    def test_log_level_verbose(self):
        "Set log level ``logging.INFO``."
        setup_logging(
            verbosity="INFO",
        )

    def test_log_level_warning(self):
        "Set log level ``logging.WARNING``."
        setup_logging(
            verbosity="WARNING",
        )

    def test_log_level_error(self):
        "Set log level ``logging.ERROR``."
        setup_logging(
            verbosity="ERROR",
        )

    def test_log_level_critical(self):
        "Set log level ``logging.CRITICAL``."
        setup_logging(
            verbosity="CRITICAL",
        )


class TestMain:
    """Test ``main()`` function."""

    def test_without_args(self):
        """Call without args."""
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 2

    def test_with_args(self, monkeypatch, tmpdir):
        """Call with args."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                str(TEST_FILE), str(TEST_FILE),
            ]
        )
        with pytest.raises(SystemExit) as exc:
            assert main() is None
        assert exc.value.code == 0

    def test_keyboard_interrupt(self, monkeypatch):

        class RaiseValueError:
            def __init__(self, *args, **kwargs):
                raise ValueError

        monkeypatch.setattr(
            'builtins.KeyboardInterrupt',
            ValueError,
        )
        monkeypatch.setattr(
            'htsinfer.cli.parse_args',
            RaiseValueError,
        )
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code >= 128


def test_main_as_script():
    """Run as script."""
    main_file = PACKAGE_DIR / "cli.py"
    fl = os.path.join(os.path.dirname(__file__), main_file)
    spec = importlib.util.spec_from_file_location('__main__', fl)
    module = importlib.util.module_from_spec(spec)
    with pytest.raises(SystemExit) as exc:
        spec.loader.exec_module(module)  # type: ignore
    assert exc.value.code == 2
