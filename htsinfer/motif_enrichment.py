"""Calculate motif enrichment from sample data."""


from scipy.stats.distributions import binom
from typing import (Dict, List)


def motif_enrichment(
    foreground: Dict[str, int],
    background: Dict[str, int]
) -> Dict:
    """Calculates enrichment and p-values of motifs with similar lengths.

    Args:
        foreground: dictionary of motifs with counts
        background: dictionary of motifs with counts

    Returns:
        Dictionary including motifs, enrichment score and p-value

    Example:
        foreground = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 2,
                      "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
        background = {"UGAUUC": 3, "UAAACC": 5, "AAGCCUUAU": 1,
                      "AGUUCUA": 1, "UUUCCCG": 5, "UUGGAA": 7}
        >>> motif_enrichment(foreground, background)
            {'UGAUUC': [2.4545454545454546, 0.0037400679729227004],
            'UAAACC': [1.0909090909090908, 0.28899727345086257],
            'AAGUUACCU': [3.2, 0.0009765625],
            'AAGCCUU': [1.8, 0.09059653269912638],
            'AGUUCUA': [0.9, 0.3899216801618356],
            'UUUCCCG': [0.9, 0.5592643397856017],
            'AAGCCUUAU': [0.26666666666666666, 0.984375],
            'UUGGAA': [0.20454545454545456, 0.9847514510038162]}
    """

    # Foreground dictionary
    pseudocode_motifs_foreground: Dict = dict()
    new_dict_foreground: Dict = dict()
    new_dict_foreground_final: Dict = dict()
    length_list_f: List = list()

    # Background dictionary
    pseudocode_motifs_background: Dict = dict()
    new_dict_background: Dict = dict()
    length_list_b: List = list()

    # Create same motifs in foreground and background dictionaries
    # Add foreground motifs to pseudocode dictionary
    for i in foreground:

        if i in background:
            pass
        else:
            pseudocode_motifs_foreground[i] = foreground[i]
            # Increase counts by amount of motfs added
            for b in background:
                if len(b) == len(i):
                    background[b] += 1
            for b in foreground:
                if len(b) == len(i):
                    foreground[b] += 1
            pseudocode_motifs_foreground[i] = 1

    # Add background motifs to pseudocode dictionary
    for i in background:

        if i in foreground:
            pass
        else:
            pseudocode_motifs_background[i] = background[i]
            # Increase counts by amount of motifs added
            for f in foreground:
                if len(f) == len(i):
                    foreground[f] += 1
            for f in background:
                if len(f) == len(i):
                    background[f] += 1
            pseudocode_motifs_background[i] = 1

    # Merge pseudocode and background
    for i in pseudocode_motifs_foreground:
        background[i] = pseudocode_motifs_foreground[i]

    # Merge pseudocode and foreground
    for i in pseudocode_motifs_background:
        foreground[i] = pseudocode_motifs_background[i]

    # Foreground dictionary
    # Calulates sum of motifs of same length
    for i in list(foreground.keys()):

        if len(i) in length_list_f:
            new_dict_foreground[len(i)] += foreground[i]
        else:
            new_dict_foreground[len(i)] = foreground[i]
            length_list_f.append(len(i))

    # Background dictionary
    # Calulates sum of motifs of same length
    for i in list(background.keys()):

        if len(i) in length_list_b:
            new_dict_background[len(i)] += background[i]
        else:
            new_dict_background[len(i)] = background[i]
            length_list_b.append(len(i))

    # Main
    # Calculates enrichment
    enrichment_dict: Dict = dict()

    for i in foreground:

        enrichment_dict[i] = ((foreground[i] * new_dict_background[len(i)]) /
                              (background[i] * new_dict_foreground[len(i)]))

    # Calculate p-value for foreground motifs
    # Create dictionary containing enrichment and p-values
    for i in foreground:
        r_bin = foreground[i]
        n_bin = new_dict_foreground[len(i)]
        p_bin = background[i] / new_dict_background[len(i)]
        new_dict_foreground_final[i] = [enrichment_dict[i],
                                        (1 - binom.cdf(r_bin, n_bin, p_bin))]

    return new_dict_foreground_final
