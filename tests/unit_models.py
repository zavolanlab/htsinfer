"""
Unit tests for '.models'.
"""
from myproj.models import (Defaults, Parameters)


def test_Parameters_init():
    res = Parameters()
    assert type(res) is Parameters


def test_Parameters_to_dict():
    res = Parameters().to_dict()
    assert type(res) is dict


def test_Defaults_init():
    res = Defaults()
    assert type(res) is Defaults


def test_Defaults_to_dict():
    res = Defaults().to_dict()
    assert type(res) is dict
