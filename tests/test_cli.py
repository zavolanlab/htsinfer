"""Unit tests for module ``cli.py``."""

import argparse
import importlib.util
import os
import sys

import pytest

from htsinfer.cli import (
    main,
    parse_args,
    setup_logging,
)
from htsinfer.models import (
    ResultsLayout,
    ResultsOrientation,
    ResultsSource,
    ResultsStats,
    ResultsType,
)
from tests.utils import (
    FILE_MATE_1,
    PACKAGE_DIR,
    RaiseError,
)


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

    def test_too_few_positional_args(self):
        """Call without positional args."""
        with pytest.raises(SystemExit) as exc:
            parse_args()
        assert exc.value.code == 2

    def test_one_positional_arg(self, monkeypatch):
        """Call with one positional args"""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                str(FILE_MATE_1),
            ]
        )
        ret_val = parse_args()
        assert isinstance(ret_val, argparse.Namespace)

    def test_two_positional_args(self, monkeypatch):
        """Call with two positional args."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                str(FILE_MATE_1), str(FILE_MATE_1),
            ]
        )
        ret_val = parse_args()
        assert isinstance(ret_val, argparse.Namespace)

    def test_too_many_positional_args(self, monkeypatch):
        """Call with too many positional args."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                str(FILE_MATE_1), str(FILE_MATE_1), str(FILE_MATE_1),
            ]
        )
        with pytest.raises(SystemExit) as exc:
            parse_args()
        assert exc.value.code == 2


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
                str(FILE_MATE_1), str(FILE_MATE_1),
                '--output-directory', str(tmpdir),
                '--temporary-directory', str(tmpdir),
            ]
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_library_source',
            ResultsSource,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_library_stats',
            ResultsStats,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_library_type',
            ResultsType,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_read_layout',
            ResultsLayout,
        )
        monkeypatch.setattr(
            'htsinfer.htsinfer.HtsInfer.get_read_orientation',
            ResultsOrientation,
        )
        with pytest.raises(SystemExit) as exc:
            assert main() is None
        assert exc.value.code == 0

    def test_with_too_many_args(self, monkeypatch):
        """Call with too many positional args."""
        monkeypatch.setattr(
            sys, 'argv', [
                'htsinfer',
                str(FILE_MATE_1), str(FILE_MATE_1), str(FILE_MATE_1),
            ]
        )
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 2

    def test_keyboard_interrupt(self, monkeypatch):
        """Test keyboard interrupt."""

        monkeypatch.setattr(
            'builtins.KeyboardInterrupt',
            ValueError,
        )
        monkeypatch.setattr(
            'htsinfer.cli.parse_args',
            RaiseError(exc=ValueError),
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
