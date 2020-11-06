"""Calculate motif enrichment from sample data."""

import numpy as np
from scipy.stats.distributions import binom
from typing import (Dict, List)

FOREGROUND_DICT = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 1,
                   "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
BACKGROUND_DICT = {"UGAUUC": 3, "UAAACC": 5, "AAGCCUUAU": 1,
                   "AGUUCUA": 1, "UUUCCCG": 5, "UUGGAA": 7}


def motif_enrichment(
    foreground: Dict[str, int],
    background: Dict[str, int]) -> Dict: 
    """Calculates enrichment and p-values of motifs with similar lengths.

    Args:
        foreground: dictionary of motifs with counts
        background: dictionary of motifs with counts

    Returns:
        Dictionary including motifs, sum of counts of motifs with
        similar lenght, enrichment score and p-value

    Example:
        foreground = {"TGATTC": 5, "TAAACC": 3, "AAGTTACCT": 1,
        "AAGCCTT": 1, "AGTTCTA": 1, "TTTCCCG": 5}
        background = {"TGATTC": 3, "TAAACC": 5, "AAGCCTTAT": 1,
        "AGTTCTA": 1, "TTTCCCG": 5, "TTGGAA": 7}
        >>> motif_enrichment(foreground, background)
            {'TGATTC': [0.7333333333333334, 0.01646090534979423],
            'TAAACC': [0.7333333333333334, 0.21947873799725662],
            'AAGTTACCT': [1.375, 0.26337448559670795],
            'AAGCCTT': [1.6041666666666667, 0.26337448559670795],
            'AGTTCTA': [1.6041666666666667, 0.26337448559670795],
            'TTTCCCG': [1.6041666666666667, 0.01646090534979423]}
    """
    # Foreground dictionary
    new_dict_foreground: Dict = dict()
    new_dict_foreground_prob: Dict = dict()
    new_dict_foreground_final: Dict = dict()
    length_list_f: List = list()

    # Background dictionary
    new_dict_background: Dict = dict()
    new_dict_background_prob: Dict = dict()
    length_list_b: List = list()

    # Foreground dictionary
    # Calulates sum of motifs of same length
    for i in list(foreground.keys()):

        if len(i) in length_list_f:
            new_dict_foreground[len(i)] += foreground[i]
        else:
            new_dict_foreground[len(i)] = foreground[i]
            length_list_f.append(len(i))

    # Calculates probability of occurrence of motifs
    for i in new_dict_foreground:
        prob_of_occurrence_f = new_dict_foreground[i] / sum(foreground.values())
        new_dict_foreground_prob[i] = prob_of_occurrence_f

    # Background dictionary
    # Calulates sum of motifs of same length
    for i in list(background.keys()):

        if len(i) in length_list_b:
            new_dict_background[len(i)] += background[i]
        else:
            new_dict_background[len(i)] = background[i]
            length_list_b.append(len(i))

    # Calculates probability of occurrence of motifs
    for i in new_dict_background:
        prob_of_occurrence_f = new_dict_background[i] / sum(background.values())
        new_dict_background_prob[i] = prob_of_occurrence_f

    # Main
    # Calculates enrichment
    enrichment_dict: Dict = dict()

    for i in new_dict_foreground_prob:
        enrichment_dict[i] = new_dict_foreground_prob[i] / new_dict_background_prob[i]

    # Calculate p-value for foreground motifs
    # Create dictionary containing enrichment and p-values
    n_binom: int = len(foreground)
    p_binom: int = np.average(list(new_dict_background_prob.values()))

    for i in foreground:
        r_binom = foreground[i]
        new_dict_foreground_final[i] = [enrichment_dict[len(i)],
                                        (1 - binom.cdf(r_binom, n_binom, p_binom))]

    return new_dict_foreground_final


motif_enrichment(FOREGROUND_DICT, BACKGROUND_DICT)
