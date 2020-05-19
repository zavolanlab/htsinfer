#!/usr/bin/env python
import importlib.util
import os

import pytest

# test parameters
CLI_FILE = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "myproj",
    "cli.py",
)
MAIN_FILE = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    "myproj",
    "__main__.py",
)
USER_INPUT = "user_input"


def test_cli(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT,
    )
    with pytest.raises(SystemExit) as error:
        fl = os.path.join(os.path.dirname(__file__), CLI_FILE)
        spec = importlib.util.spec_from_file_location('__main__', fl)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    assert error.value.code == 0


def test_main(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT,
    )
    with pytest.raises(SystemExit) as error:
        fl = os.path.join(os.path.dirname(__file__), MAIN_FILE)
        spec = importlib.util.spec_from_file_location('__main__', fl)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    assert error.value.code == 0
