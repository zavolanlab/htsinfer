"""
Unit tests for '.params'.
"""
from string import (ascii_lowercase, ascii_uppercase, digits)

import pytest

from myproj.params import GetParams
from myproj.models import (Defaults, Parameters)

# Test parameters
DEFAULTS = Defaults().to_dict()
PARAMS = Parameters().to_dict()
UNSANITIZED = "@#^;?/My_!!!$pRO-JE*CT  123&&    &"
LIST = [1, 2, 3]
SLUG_DEFAULT = "my_project_123"
SLUG_UPPERCASE_REMOVED = "y_p_123"
SLUG_KEEP_UPPERCASE = "My_pROJECT_123"
SLUG_WHITESPACE_REMOVED = "my_project123"
DEFAULTS_DICT_CORRUPT = {}
DEFAULTS_DICT_GENERIC = {
    "value": "value",
    "description": "description",
}
DEFAULTS_DICT_CHOICES = {
    "value": "value",
    "choices": ['a', 'b', 'c'],
    "multiple": False,
    "alternative": 'none',
    "description": "description",
}
DEFAULTS_DICT_CHOICES_MULTI = {
    "value": "value",
    "choices": ['a', 'b', 'c'],
    "multiple": True,
    "alternative": 'none',
    "description": "description",
}
DEFAULTS_DICT_CHOICES_NO_ALT = {
    "value": "value",
    "choices": ['a', 'b', 'c'],
    "multiple": False,
    "alternative": None,
    "description": "description",
}
USER_INPUT_GENERIC = "user_input"
USER_INPUT_CHOICES = "a"
USER_INPUT_CHOICES_INVALID = "not_a_choice"
RETRIES = 2
SEP_MULTI = [',', '|', '&']
CHOICES = "a"
CHOICES_SPLIT = ["a"]
CHOICES_MULTI = " a & b, c &  d | e "
CHOICES_MULTI_SPLIT = ["a & b", "c &  d | e"]
CHOICES_MULTI_SPLIT_MULTI_SEP = ["a", "b", "c", "d", "e"]


# __init__()
def test_init_no_args():
    with pytest.raises(TypeError):
        GetParams()


def test_init_args_all_params():
    p = GetParams(
        defaults=DEFAULTS,
        params=PARAMS,
    )
    assert p.params == PARAMS


# get_params()
def test_get_params_no_org(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    p = GetParams(
        defaults=DEFAULTS,
        params=PARAMS,
    )
    del p.params['org']
    assert p.get_params(retries=1) is None


def test_get_params_no_user(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    p = GetParams(
        defaults=DEFAULTS,
        params=PARAMS,
    )
    del p.params['user']
    assert p.get_params(retries=1) is None


def test_get_params_popen_exception(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    monkeypatch.setattr(
        'os.popen',
        lambda description: (_ for _ in ()).throw(Exception),
    )
    p = GetParams(
        defaults=DEFAULTS,
        params=PARAMS,
    )
    del p.params['user']
    assert p.get_params(retries=1) is None


def test_get_params_no_project(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    p = GetParams(
        defaults=DEFAULTS,
        params=PARAMS,
    )
    del p.params['project']
    assert p.get_params(retries=1) is None


def test_get_params_no_soft(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    p = GetParams(
        defaults=DEFAULTS,
        params=PARAMS,
    )
    del p.params['soft']
    assert p.get_params(retries=1) is None


# query_user()
def test_query_user_no_args(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    with pytest.raises(TypeError):
        GetParams.query_user()


def test_query_user_wrong_type_d():
    with pytest.raises(TypeError):
        GetParams.query_user(
            d=LIST,
        )


def test_query_user_wrong_type_retries():
    with pytest.raises(TypeError):
        GetParams.query_user(
            d=DEFAULTS_DICT_CORRUPT,
            retries=LIST,
        )


def test_query_user_defaults_dict_corrupt(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    with pytest.raises(TypeError):
        GetParams.query_user(
            d=DEFAULTS_DICT_CORRUPT,
        )


def test_query_user_defaults_dict_generic(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_GENERIC,
    )
    assert GetParams.query_user(
        d=DEFAULTS_DICT_GENERIC,
    ) == USER_INPUT_GENERIC


def test_query_user_defaults_dict_generic_no_input(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: "",
    )
    assert GetParams.query_user(
        d=DEFAULTS_DICT_GENERIC,
    ) == DEFAULTS_DICT_GENERIC['value']


def test_query_user_defaults_dict_choices(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_CHOICES,
    )
    assert GetParams.query_user(
        d=DEFAULTS_DICT_CHOICES
    ) == USER_INPUT_CHOICES


def test_query_user_defaults_dict_choices_multi(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_CHOICES,
    )
    assert GetParams.query_user(
        d=DEFAULTS_DICT_CHOICES_MULTI
    ) == USER_INPUT_CHOICES


def test_query_user_defaults_dict_choices_no_alt(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_CHOICES,
    )
    assert GetParams.query_user(
        d=DEFAULTS_DICT_CHOICES_NO_ALT
    ) == USER_INPUT_CHOICES


def test_query_user_defaults_dict_choices_wrong_input(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_CHOICES_INVALID,
    )
    with pytest.raises(ValueError):
        GetParams.query_user(
            d=DEFAULTS_DICT_CHOICES,
            retries=RETRIES,
        ) == USER_INPUT_CHOICES


def test_query_user_defaults_dict_choices_multi_wrong_input(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: USER_INPUT_CHOICES_INVALID,
    )
    with pytest.raises(ValueError):
        GetParams.query_user(
            d=DEFAULTS_DICT_CHOICES_MULTI,
            retries=RETRIES,
        ) == USER_INPUT_CHOICES


def test_query_user_raise_keyboard_interrupt(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda description: (_ for _ in ()).throw(KeyboardInterrupt),
    )
    with pytest.raises(KeyboardInterrupt):
        GetParams.query_user(
            d=DEFAULTS_DICT_GENERIC,
        )


# slugify()
def test_slugify_no_args():
    with pytest.raises(TypeError):
        GetParams.slugify()


def test_slugify_okay():
    assert GetParams.slugify(
        s=UNSANITIZED,
    ) == SLUG_DEFAULT


def test_slugify_wrong_type_s():
    with pytest.raises(TypeError):
        GetParams.slugify(
            s=LIST,
        )


def test_slugify_wrong_type_allowed_start():
    with pytest.raises(TypeError):
        GetParams.slugify(
            s=UNSANITIZED,
            allowed_start=LIST,
        )


def test_slugify_wrong_type_allowed_end():
    with pytest.raises(TypeError):
        GetParams.slugify(
            s=UNSANITIZED,
            allowed_end=LIST,
        )


def test_slugify_wrong_type_allowed_rest():
    with pytest.raises(TypeError):
        GetParams.slugify(
            s=UNSANITIZED,
            allowed_rest=LIST,
        )


def test_slugify_wrong_type_lower():
    with pytest.raises(TypeError):
        GetParams.slugify(
            s=UNSANITIZED,
            lower=LIST,
        )


def test_slugify_wrong_type_whitespace_replace():
    with pytest.raises(TypeError):
        GetParams.slugify(
            s=UNSANITIZED,
            whitespace_replace=LIST,
        )


def test_slugify_no_whitespace():
    assert GetParams.slugify(
        s=UNSANITIZED,
        whitespace_replace="",
    ) == SLUG_WHITESPACE_REMOVED


def test_slugify_remove_uppercase():
    assert GetParams.slugify(
        s=UNSANITIZED,
        lower=False,
    ) == SLUG_UPPERCASE_REMOVED


def test_slugify_keep_uppercase():
    assert GetParams.slugify(
        s=UNSANITIZED,
        allowed_start=ascii_uppercase + ascii_lowercase,
        allowed_end=ascii_uppercase + ascii_lowercase + digits,
        allowed_rest=ascii_uppercase + ascii_lowercase + digits + "_",
        lower=False,
    ) == SLUG_KEEP_UPPERCASE


def test_slugify_no_chars_allowed_at_start():
    assert GetParams.slugify(
        s=UNSANITIZED,
        allowed_start="",
    ) == ""


def test_slugify_no_chars_allowed_at_end():
    assert GetParams.slugify(
        s=UNSANITIZED,
        allowed_end="",
    ) == ""


# split_choices()
def test_split_choices_no_args():
    with pytest.raises(TypeError):
        GetParams.split_choices()


def test_split_choices_wrong_type_s():
    with pytest.raises(TypeError):
        GetParams.split_choices(
            s=LIST,
        )


def test_split_choices_wrong_type_sep_1():
    with pytest.raises(TypeError):
        GetParams.split_choices(
            s=UNSANITIZED,
            sep=UNSANITIZED,
        )


def test_split_choices_wrong_type_sep_2():
    with pytest.raises(TypeError):
        GetParams.split_choices(
            s=UNSANITIZED,
            sep=LIST,
        )


def test_split_choices_one_item():
    assert GetParams.split_choices(
        s=CHOICES,
    ) == CHOICES_SPLIT


def test_split_choices_multi_items():
    assert GetParams.split_choices(
        s=CHOICES_MULTI,
    ) == CHOICES_MULTI_SPLIT


def test_split_choices_multi_items_multi_sep():
    assert GetParams.split_choices(
        s=CHOICES_MULTI,
        sep=SEP_MULTI,
    ) == CHOICES_MULTI_SPLIT_MULTI_SEP
