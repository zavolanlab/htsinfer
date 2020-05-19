"""
Unit tests for '.config'.
"""
import os

import pytest
from yaml.parser import ParserError
from yaml.representer import RepresenterError

from myproj.config import ConfigParser
from myproj.models import Parameters

# Test parameters
FILE_OK = os.path.join(
    os.path.dirname(__file__),
    "files",
    "yaml",
)
FILE_UNAVAILABLE = "xyz/zyx/123"
FILE_NOT_YAML = __file__
FILE_EMPTY = os.path.join(
    os.path.dirname(__file__),
    "files",
    "empty",
)
FILE_TXT = os.path.join(
    os.path.dirname(__file__),
    "files",
    "txt",
)
FILE_OUT = os.path.join(
    os.path.dirname(__file__),
    "files",
    "conf_out",
)
STRING = "SOME HEADER"
KWARGS = {
    "STRING": STRING,
    "INTEGER": 123,
    "DICT": {"abc": 1, "cde": 2, "efg": 3},
    "DICT_EMPTY": {},
    "LIST": [1, 2, 3],
    "LIST_EMPTY": [],
}
KEY_1 = "a"
KEY_2 = "a1"
KEY_3 = "a2"
KEY_4 = "b"
KEY_5 = "c"
INT = 1
LIST = [1, 2, 3]
OBJECT = {"OBJECT": ConfigParser}
DICT_1 = {KEY_1: {KEY_2: 2, KEY_3: 3}}
DICT_2 = {KEY_1: {KEY_2: 5}, KEY_4: 6}
QUERY = {KEY_1: {KEY_2: INT, KEY_3: {}}, KEY_4: INT, KEY_5: KEY_1}
QUERY_FALSE = {KEY_1: INT, KEY_4: INT, KEY_5: KEY_1}
REF = {KEY_1: {KEY_2: INT}, KEY_4: [], KEY_5: {}}


# __init__()
def test_init_no_args():
    res = ConfigParser()
    assert res.values == {}


def test_init_single_config():
    res = ConfigParser(FILE_OK)
    assert type(res.values) is dict


def test_init_config_unavailable():
    with pytest.raises(FileNotFoundError):
        ConfigParser(FILE_UNAVAILABLE)


def test_init_config_invalid():
    with pytest.raises(ParserError):
        ConfigParser(FILE_NOT_YAML)


def test_init_multi_config():
    res = ConfigParser(FILE_OK, FILE_OK)
    assert type(res.values) is dict


def test_init_single_config_log():
    res = ConfigParser(FILE_OK, log=True)
    assert type(res.values) is dict


def test_init_empty_config_log():
    res = ConfigParser(FILE_EMPTY, log=True)
    assert res.values == {}


# log_yaml()
def test_log_yaml_no_args():
    assert ConfigParser.log_yaml() is None


def test_log_yaml_header():
    assert ConfigParser.log_yaml(header=STRING) is None


def test_log_yaml_kwargs():
    assert ConfigParser.log_yaml(kwargs=KWARGS) is None


def test_log_yaml_kwargs_with_error():
    with pytest.raises(TypeError):
        ConfigParser.log_yaml(**KWARGS, level=OBJECT)


def test_log_yaml_header_and_kwargs():
    assert ConfigParser.log_yaml(header=STRING, **KWARGS) is None


# read_config_files()
def test_read_config_files_no_args():
    res = ConfigParser.read_config_files()
    assert res == {}


def test_read_config_files_single_config():
    res = ConfigParser.read_config_files(FILE_OK)
    assert type(res) is dict


def test_read_config_files_multi_config():
    res = ConfigParser.read_config_files(FILE_OK, FILE_OK)
    assert type(res) is dict


def test_read_config_files_single_config_unavailable():
    with pytest.raises(FileNotFoundError):
        ConfigParser.read_config_files(FILE_UNAVAILABLE)


def test_read_config_files_multi_config_partly_unavailable():
    with pytest.raises(FileNotFoundError):
        ConfigParser.read_config_files(FILE_OK, FILE_UNAVAILABLE)


def test_read_config_files_single_config_invalid():
    with pytest.raises(ParserError):
        ConfigParser.read_config_files(FILE_NOT_YAML)


def test_read_config_files_multi_config_partly_invalid():
    with pytest.raises(ParserError):
        ConfigParser.read_config_files(FILE_OK, FILE_NOT_YAML)


# recursive_dict_update()
def test_recursive_dict_update_correct_inputs():
    d = ConfigParser.recursive_dict_update(
        original=DICT_1,
        update=DICT_2,
    )
    assert d[KEY_1][KEY_2] == DICT_2[KEY_1][KEY_2]
    assert KEY_3 in d[KEY_1]
    assert KEY_4 in d


def test_recursive_dict_update_arg1_list():
    with pytest.raises(TypeError):
        ConfigParser.recursive_dict_update(
            original=LIST,
            update=DICT_1,
        )


def test_recursive_dict_update_arg2_list():
    with pytest.raises(TypeError):
        ConfigParser.recursive_dict_update(
            original=DICT_1,
            update=LIST,
        )


def test_recursive_dict_update_arg1_int():
    with pytest.raises(TypeError):
        ConfigParser.recursive_dict_update(
            original=INT,
            update=DICT_1,
        )


def test_recursive_dict_update_arg2_int():
    with pytest.raises(TypeError):
        ConfigParser.recursive_dict_update(
            original=DICT_1,
            update=INT,
        )


def test_recursive_dict_update_arg1_str():
    with pytest.raises(TypeError):
        ConfigParser.recursive_dict_update(
            original=KEY_1,
            update=DICT_1,
        )


def test_recursive_dict_update_arg2_str():
    with pytest.raises(TypeError):
        ConfigParser.recursive_dict_update(
            original=DICT_1,
            update=KEY_1,
        )


# same_keys()
def test_same_keys_correct_inputs():
    assert ConfigParser.same_keys(
        query=QUERY,
        ref=REF,
    ) is True


def test_same_keys_correct_inputs_two_way():
    assert ConfigParser.same_keys(
        query=QUERY,
        ref=QUERY,
        two_way=True,
    ) is True


def test_same_keys_conflicting_inputs():
    assert ConfigParser.same_keys(
        query=QUERY_FALSE,
        ref=REF,
    ) is False


def test_same_keys_conflicting_inputs_two_way():
    assert ConfigParser.same_keys(
        query=QUERY,
        ref=REF,
        two_way=True,
    ) is False


def test_same_keys_wrong_types_ref():
    with pytest.raises(TypeError):
        ConfigParser.same_keys(
            query=QUERY,
            ref=INT,
        )


def test_same_keys_wrong_types_query_list():
    assert ConfigParser.same_keys(
            query=LIST,
            ref=REF,
    ) is False


def test_same_keys_wrong_types_query_int():
    assert ConfigParser.same_keys(
            query=INT,
            ref=REF,
    ) is False


def test_same_keys_wrong_types_query_str():
    assert ConfigParser.same_keys(
            query=KEY_1,
            ref=REF,
    ) is False


def test_same_keys_wrong_types_query_none():
    assert ConfigParser.same_keys(
            query=None,
            ref=REF,
    ) is False


def test_same_keys_wrong_types_query_class():
    assert ConfigParser.same_keys(
            query=ConfigParser,
            ref=REF,
    ) is False


# dict_to_yaml()
def test_dict_to_yaml_file_ok():
    params = Parameters().to_dict()
    ret = ConfigParser.dict_to_yaml(
        d=params,
        yaml_file=FILE_OUT,
    )
    assert ret is None


def test_dict_to_yaml_wrong_type_d():
    with pytest.raises(TypeError):
        ConfigParser.dict_to_yaml(
            d=LIST,
            yaml_file=FILE_OUT,
        )


def test_dict_to_yaml_wrong_type_yaml_file():
    params = Parameters().to_dict()
    with pytest.raises(TypeError):
        ConfigParser.dict_to_yaml(
            d=params,
            yaml_file=LIST,
        )


def test_dict_to_yaml_file_unavailable():
    params = Parameters().to_dict()
    ret = ConfigParser.dict_to_yaml(
        d=params,
        yaml_file=FILE_UNAVAILABLE,
    )
    assert ret is None


def test_dict_to_yaml_invalid_object():
    params = {Parameters(): Parameters()}
    with pytest.raises(RepresenterError):
        ConfigParser.dict_to_yaml(
            d=params,
            yaml_file=FILE_OUT,
        )


# yaml_to_dict()
def test_yaml_to_dict_file_ok():
    d = ConfigParser.yaml_to_dict(yaml_file=FILE_OK)
    assert type(d) is dict
    assert bool(d) is True


def test_yaml_to_dict_file_not_found():
    with pytest.raises(FileNotFoundError):
        ConfigParser.yaml_to_dict(yaml_file=FILE_UNAVAILABLE)


def test_yaml_to_dict_file_not_yaml():
    with pytest.raises(ParserError):
        ConfigParser.yaml_to_dict(yaml_file=FILE_NOT_YAML)


def test_yaml_to_dict_file_txt():
    with pytest.raises(TypeError):
        ConfigParser.yaml_to_dict(yaml_file=FILE_TXT)


def test_yaml_to_dict_file_empty():
    assert ConfigParser.yaml_to_dict(yaml_file=FILE_EMPTY) == {}
