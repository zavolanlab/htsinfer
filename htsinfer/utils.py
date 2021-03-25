"""Validator function for infer_organism and infer_adapter"""

from pandas import DataFrame  # type: ignore


def minmatch_factor_validator(
    count_df: DataFrame,
    column_index: int,
    min_match: float = 5,
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
            considered as the resulting adapter/organism.

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
