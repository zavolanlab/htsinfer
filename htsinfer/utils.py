"""Utilities used across multiple HTSinfer modules."""

from pandas import DataFrame  # type: ignore


def results_validator(
    data: DataFrame,
    column_index: int,
    min_value: float = 5,
    min_ratio: float = 2
) -> bool:
    """Validates whether (1) the first value of a data frame column is equal to
    or higher than a specified minimum value, and (2) that the ratio of the
    first and second values of the column is higher than a specified minimum
    ratio.

    Args:
        data: Pandas data frame containing values to compare in a single
            column.
        column_index: Index of `data` column containing values of interest.
        min_value: Minimum value required in first row of `column_index` for
            validation to pass.
        min_ratio: Minimum ratio of first and second rows of `column_index`
            required for validation to pass.

    Returns:
        Whether data frame `data` satisfies the `min_value` and `min_ratio`
        constraints for value in column `column_index`.
    """
    if data.iloc[0][column_index] < min_value:
        return False
    if data.iloc[1][column_index] != 0:
        ratio = data.iloc[0][column_index]/data.iloc[1][column_index]
        if ratio < min_ratio:
            return False
    return True
