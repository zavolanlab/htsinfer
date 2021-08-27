"""Unit tests for module ``utils.py``."""

import pytest

from htsinfer.utils import (
    convert_dict_to_df,
    validate_top_score,
)
from tests.utils import DICT_DF


class TestConvertDictToDf:
    """Test ``convert_dict_to_df()`` function."""

    COLUMN_HEADERS = ('COL_1', 'COL_2')

    def test_error_too_few_headers(self):
        """Test invalid (too few) column headers input."""
        with pytest.raises(ValueError):
            convert_dict_to_df(
                dic=DICT_DF,
                col_headers=('COL_1')  # type: ignore
            )

    def test_error_too_many_headers(self):
        """Test invalid (too many) column headers input."""
        with pytest.raises(ValueError):
            convert_dict_to_df(
                dic=DICT_DF,
                col_headers=('COL_1', 'COL_2', 'COL_3')  # type: ignore
            )

    def test_valid(self):
        """Test valid conversion."""
        df = convert_dict_to_df(
            dic=DICT_DF,
        )
        assert df.shape == (3, 2)

    def test_valid_headers(self):
        """Test valid conversion with headers."""
        df = convert_dict_to_df(
            dic=DICT_DF,
            col_headers=self.COLUMN_HEADERS,
        )
        assert df.shape == (3, 2)
        assert tuple(df.columns.values) == self.COLUMN_HEADERS

    def test_valid_sort(self):
        """Test valid conversion with sorting."""
        df = convert_dict_to_df(
            dic=DICT_DF,
            sort=True,
        )
        assert df.shape == (3, 2)
        assert df.iloc[:, 0].to_list() == sorted(DICT_DF.keys())
        assert df.at[0, 0] == 'key_1'

    def test_valid_sort_rev(self):
        """Test valid conversion with reverse sorting."""
        df = convert_dict_to_df(
            dic=DICT_DF,
            sort=True,
            sort_ascending=False,
        )
        assert df.shape == (3, 2)
        assert df.iloc[:, 0].to_list() == sorted(DICT_DF.keys(), reverse=True)
        assert df.at[0, 0] == 'key_3'

    def test_valid_sort_by_index(self):
        """Test valid conversion with sorting by non-default column."""
        df = convert_dict_to_df(
            dic=DICT_DF,
            col_headers=self.COLUMN_HEADERS,
            sort=True,
            sort_by=1,
        )
        assert df.shape == (3, 2)
        assert df.iloc[:, 1].to_list() == sorted(DICT_DF.values())
        assert df.at[0, 'COL_2'] == 1


class TestValidateTopScore:
    """Test ``validate_top_score()`` function."""

    def test_error_input(self):
        """Test invalid input."""
        with pytest.raises(ValueError):
            validate_top_score(
                vector=[100, 'a', 0, 0],  # type: ignore
            )

    def test_valid(self):
        """Test input fulfilling all validation criteria."""
        res = validate_top_score(
            vector=[100, 10, 1],
        )
        assert res is True

    def test_valid_div_zero(self):
        """Test input fulfilling all validation criteria, with second highest
        value being zero.
        """
        res = validate_top_score(
            vector=[100, 0, 0],
        )
        assert res is True

    def test_valid_sort(self):
        """Test input fulfilling all validation criteria after sorting."""
        res = validate_top_score(
            vector=[1, 10, 100],
            rev_sorted=False,
        )
        assert res is True

    def test_invalid_div_zero(self):
        """Test input failing validation due to the second highest value being
        zero.
        """
        res = validate_top_score(
            vector=[100, 0, 0],
            accept_zero=False,
        )
        assert res is False

    def test_invalid_unsorted(self):
        """Test input failing validation due to sorting."""
        res = validate_top_score(
            vector=[1, 10, 100],
            rev_sorted=True,
        )
        assert res is False

    def test_invalid_value_too_small(self):
        """Test input failing validation due to the top value not exceeding
        the minimum value.
        """
        res = validate_top_score(
            vector=[4, 1, 0],
            min_value=5,
        )
        assert res is False

    def test_invalid_ratio_too_small(self):
        """Test input failing validation due to the top value not being
        sufficiently bigger than the second highest value.
        """
        res = validate_top_score(
            vector=[4, 3, 0],
            min_ratio=2,
        )
        assert res is False

    def test_invalid_no_items(self):
        """Test input failing because the input list is empty."""
        res = validate_top_score(
            vector=[],
        )
        assert res is False

    def test_invalid_one_item(self):
        """
        Test input failing because the input list contains only a single
        item.
        """
        res = validate_top_score(
            vector=[1000],
        )
        assert res is False
