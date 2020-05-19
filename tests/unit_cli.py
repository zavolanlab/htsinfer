"""
Unit tests for module '.cli'.
"""
import logging
import os

import pytest

from myproj.cli import (main, parse_cli_args, setup_logging)

# test parameters
CONFIG_FILE = os.path.join(
    os.path.dirname(__file__),
    "files",
    "yaml"
)
DEFAULTS = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "config",
    "defaults.yaml",
)
PARAMS = os.path.join(
    os.path.dirname(__file__),
    "files",
    "params",
)
PARAMS_EXTRA = os.path.join(
    os.path.dirname(__file__),
    "files",
    "params_extra",
)
HELP_OPTION = "help"
SWITCH = "debug"
SWITCH_DEFAULT = False
INVALID_OPTION = "--zyxw123"
FILE_OPTION = "defaults"
MULTI_FILE_OPTION = "config"
VALID_FILE = __file__
INVALID_FILE = "/xyz/zyx/125"
LOGGER = logging.getLogger()
DEFAULT_LOG_LEVEL = logging.WARNING
VERBOSE_LOG_LEVEL = logging.INFO
DEBUG_LOG_LEVEL = logging.DEBUG
USER_INPUT = "user_input"
DEFAULT_DICT_GENERIC = {'value': "", 'description': ""}


# main()
def test_main_with_default_args(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT,
    )
    with pytest.raises(SystemExit) as e:
        main()
    assert e.value.code == 0


def test_main_with_no_config_files(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT,
    )
    with pytest.raises(SystemExit) as e:
        main(
            config_files=None,
        )
    assert e.value.code == 0


def test_main_with_explicit_args(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT,
    )
    with pytest.raises(SystemExit) as e:
        main(
            defaults_file=DEFAULTS,
            config_files=[CONFIG_FILE],
        )
    assert e.value.code == 0


def test_main_with_config_str():
    with pytest.raises(SystemExit) as e:
        main(
            config_files=CONFIG_FILE,
        )
    assert e.value.code == 1


def test_main_with_corrupt_defaults():
    with pytest.raises(SystemExit) as e:
        main(
            defaults_file=CONFIG_FILE,
        )
    assert e.value.code == 1


def test_main_with_defaults_and_complete_config_list():
    with pytest.raises(SystemExit) as e:
        main(
            defaults_file=DEFAULTS,
            config_files=[PARAMS],
        )
    assert e.value.code == 0


# parse_cli_args()
def test_help_option():
    with pytest.raises(SystemExit):
        assert parse_cli_args(["--" + HELP_OPTION])


def test_no_args():
    ret = parse_cli_args([])
    assert vars(ret)[SWITCH] is SWITCH_DEFAULT


def test_switch():
    ret = parse_cli_args(["--" + SWITCH])
    assert vars(ret)[SWITCH] is not SWITCH_DEFAULT


def test_invalid_option():
    with pytest.raises(SystemExit):
        assert parse_cli_args([INVALID_OPTION])


def test_file_returns_string():
    ret = parse_cli_args(["--" + FILE_OPTION, VALID_FILE])
    assert type(vars(ret)[FILE_OPTION]) is str


def test_multi_file_returns_list_of_strings():
    ret = parse_cli_args([
        "--" + MULTI_FILE_OPTION,
        VALID_FILE,
        "--" + MULTI_FILE_OPTION,
        VALID_FILE
    ])
    assert type(vars(ret)[MULTI_FILE_OPTION]) is list
    for arg in vars(ret)[MULTI_FILE_OPTION]:
        assert type(arg) is str


def test_action_open_invalid_file():
    with pytest.raises(SystemExit):
        assert parse_cli_args(["--" + FILE_OPTION, INVALID_FILE])


def test_no_arg():
    with pytest.raises(SystemExit):
        assert parse_cli_args(["--" + FILE_OPTION])


# setup_logging()
def test_default_log_level():
    setup_logging(LOGGER)
    assert LOGGER.level == DEFAULT_LOG_LEVEL


def test_verbose_log_level():
    setup_logging(LOGGER, verbose=True)
    assert LOGGER.level == VERBOSE_LOG_LEVEL


def test_debug_log_level():
    setup_logging(LOGGER, debug=True)
    assert LOGGER.level == DEBUG_LOG_LEVEL


def test_log_level_precedence():
    setup_logging(LOGGER, verbose=True, debug=True)
    assert LOGGER.level == DEBUG_LOG_LEVEL


def test_missing_logger():
    assert setup_logging() is None
