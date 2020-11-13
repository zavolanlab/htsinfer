"""Calculate motif enrichment from sample data."""

from typing import (Dict, List)
from scipy.stats import binom


def motif_adding_foreground(
        foreground: Dict[str, int],
        background: Dict[str, int]
) -> Dict:
    """Creates dictionary with motifs included in foreground but
        not background. Increases count of motifs with similar lengths
        in foreground by one for each motif missing in background.
        Adds pseudo-code of 1 to each motif added to background.
        Increases count of motifs with similar lengths
        in background by one for each added motif.

    Args:
        foreground: dictionary of motifs with counts
        background: dictionary of motifs with counts

    Returns:
        Dictionary with missing motifs in background.

    Example:
        foreground = {"UGAUUC": 5, "CCUUAA": 3, "AAGUUACCU": 2,
                      "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
        background = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                      "AAGCCUU": 1, "AGUUCUA": 4, "UUUAAG": 3}
        >>> motif_adding(foreground, background)
           {'CCUUAA': 1, 'UUUCCCG': 1}
    """
    # Pseudocode dictionary foreground
    pseudocode_motifs_foreground: Dict = dict()

    # Add foreground motifs to pseudocode dictionary
    for i in foreground:

        if i in background:
            pass
        else:
            pseudocode_motifs_foreground[i] = foreground[i]
            # Increase counts by amount of motfs added
            for value in background:
                if len(value) == len(i):
                    background[value] += 1
            for value in foreground:
                if len(value) == len(i):
                    foreground[value] += 1
            pseudocode_motifs_foreground[i] = 1

    return pseudocode_motifs_foreground


def motif_adding_background(
        foreground: Dict[str, int],
        background: Dict[str, int]
) -> Dict:
    """Creates dictionary with motifs included in background but
        not foreground. Increases count of motifs with similar lengths
        in background by one for each motif missing in foreground.
        Adds pseudo-code of 1 to each motif added to foreground.
        Increases count of motifs with similar lengths
        in foreground by one for each added motif.

    Args:
        foreground: dictionary of motifs with counts
        background: dictionary of motifs with counts

    Returns:
        Dictionary with motifs missing in foreground.

    Example:
        foreground = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 2,
                      "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
        background = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                      "AAGCCUU": 1, "AGUUCUA": 4, "UUUCCCG": 3}
        >>> motif_adding_background(foreground, background)
            {'UAAACC': 1, 'UUUAAG': 1}
    """
    # Pseudocode dictionary background
    pseudocode_motifs_background: Dict = dict()

    # Add background motifs to pseudocode dictionary
    for i in background:

        if i in foreground:
            pass
        else:
            pseudocode_motifs_background[i] = background[i]
            # Increase counts by amount of motifs added
            for value in foreground:
                if len(value) == len(i):
                    foreground[value] += 1
            for value in background:
                if len(value) == len(i):
                    background[value] += 1
            pseudocode_motifs_background[i] = 1

    return pseudocode_motifs_background


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
        background = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                      "AAGCCUU": 1, "AGUUCUA": 4, "UUUCCCG": 3}
        >>> motif_enrichment(foreground, background)
            {'UGAUUC': [0.6944444444444444, 0.96190821],
             'UAAACC': [3.75, 0.005024350000000011],
             'AAGUUACCU': [1.0, 0.0],
             'AAGCCUU': [1.1428571428571428, 0.21460819244384766],
             'AGUUCUA': [0.2857142857142857, 0.9375],
             'UUUCCCG': [1.9047619047619047, 0.013209342956542969]}
    """
    # Create same motifs in both dictionaries
    pseudo_foreground = motif_adding_foreground(foreground, background)
    pseudo_background = motif_adding_background(foreground, background)

    # Merge pseudocode and dictionaries
    for i in pseudo_background:
        foreground[i] = pseudo_background[i]

    for i in pseudo_foreground:
        background[i] = pseudo_foreground[i]

    # Foreground dictionary
    new_dict_foreground: Dict = dict()
    new_dict_foreground_final: Dict = dict()
    length_list_f: List = list()

    # Background dictionary
    new_dict_background: Dict = dict()
    length_list_b: List = list()

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
        n_bin = new_dict_foreground[len(i)]
        p_bin = background[i] / new_dict_background[len(i)]
        new_dict_foreground_final[i] = [enrichment_dict[i],
                                        (1 - binom.cdf(foreground[i],
                                                       n_bin, p_bin))]

    return new_dict_foreground_final
