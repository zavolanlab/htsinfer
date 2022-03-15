"""Utilities used across multiple HTSinfer modules."""
from typing import (Dict, List, Optional, Tuple)
from pandas import DataFrame  # type: ignore


def convert_dict_to_df(
    dic: Dict,
    col_headers: Optional[Tuple[str, str]] = None,
    sort: bool = False,
    sort_by: int = 0,
    sort_ascending: bool = True,
) -> DataFrame:
    """Convert dictionary to two-column data frame.

    Args:
        dic: Dictionary to convert.
        col_headers: List of column headers. Length MUST match number of
            dictionary keys/data frame columns.
        sort: Whether the resulting data frame is supposed to be sorted.
        sort_by: Column index used for sorting. Ignored if `sort` is `False`.
        sort_ascending: Whether the data frame is supposed to be sorted in
            ascending order. Ignored if `sort` is `False`.

    Returns:
        Data frame prepared from dictionary.

    Raises:
        ValueError: Raised if number of provided column headers does not match
            the number of data frame columns.
    """
    dat = DataFrame(dic.items())
    if col_headers is not None:
        if len(col_headers) != len(dat.columns):
            raise ValueError(
                f"Number of column headers ({len(col_headers)}) does not "
                f"match number of data frame columns ({len(dat.columns)})."
            )
        dat.columns = list(col_headers)
    if sort:
        dat = dat.sort_values(
            by=dat.columns[sort_by],
            ascending=sort_ascending,
        ).reset_index(drop=True)
    return dat


def validate_top_score(
    vector: List[float],
    min_value: float = 2,
    min_ratio: float = 2,
    accept_zero: bool = True,
    rev_sorted: bool = True,
) -> bool:
    """Validates whether (1) the maximum value of a numeric list is equal to or
    higher than a specified minimum value AND (2) that the ratio of the first
    and second highest values of the list is higher than a specified minimum
    ratio.

    If the passed list/vector does NOT contain at least two items, the function
    returns `False`.

    Args:
        vector: List of numbers.
        min_value: Minimum value required in first row of `column_index` for
            validation to pass.
        min_ratio: Minimum ratio of first and second rows of `column_index`
            required for validation to pass.
        accept_zero: Whether to accept a top score (i.e., return `True`) if the
            second highest value in the provided list is zero. If not set to
            `True`, `False` is returned in these cases.
        rev_sorted: Whether the list of numbers is sorted in descencing numeric
            order.

    Returns:
        Whether data frame `data` satisfies the `min_value` and `min_ratio`
        constraints for value in column `column_index`.

    Raises:
        ValueError: Raised if one of the list items can not be interpreted as a
            number.
    """
    try:
        vector = [float(item) for item in vector]
    except ValueError as exc:
        raise ValueError(
            f"Input list '{vector}' contains non-numeric items."
        ) from exc
    if not rev_sorted:
        vector = sorted(vector, reverse=True)
    try:
        if vector[0] < min_value:
            return False
        if vector[0] / vector[1] < min_ratio:
            return False
    except IndexError:
        return False
    except ZeroDivisionError:
        if not accept_zero:
            return False
    return True


def minmatch_factor_validator(
    count_df: DataFrame,
    column_index: int,
    min_match: float = 10,
    factor: float = 2
) -> bool:
    """Validates min_match and factor for adapter/organism to be considered as
    resulting adapter/organism.

    Args:
        count_df:
            adapter: Adapters sequence count percentage.
            organism: Count info percentage for all organisms.
        column_index: Column parameter for adapter and organism dataframe.
        min_match: Minimum percentage of reads that contain a given
            adapter/organism in order for that adapter/organism to be
            considered as the resulting adapter/organism.
        factor: The minimum frequency ratio between the first and second most
            frequent adapter/organism in order for an adapter/organism to be
            considered as resulting adapter/organism.

    Returns:
        Whether it satisfies the minimum match percentage and the minimum
        frequency ratio.
    """
    if count_df.iloc[0][column_index] < min_match:
        return False
    if count_df.iloc[1][column_index] != 0:
        ratio = count_df.iloc[0][column_index]/count_df.iloc[1][column_index]
        if ratio < factor:
            return False
    return True
