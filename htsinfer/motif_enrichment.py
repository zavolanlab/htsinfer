"""Calculate motif enrichment from sample data."""

from typing import (Dict, List)
from scipy.stats import binom

def motif_enrichment(
        foreground: Dict[str, int],
        background: Dict[str, int]
) -> Dict:
    """Calculates enrichment and p-values of motifs in foreground with respect
    to background sample. 
    Calculation done separately for distinct motif lengths.

    Args:
        foreground: dictionary of motifs with counts
        background: dictionary of motifs with counts

    Returns:
        Dictionary with motifs as keys and enrichment score and p-value as values

    Example:
        foreground = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 2,
                      "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
        background = {"UGAUUC": 9, "UAAACC": 1, "AAGUUACCU": 3,
                      "AAGCCUU": 1, "AGUUCUA": 4, "UUUCCCG": 3}
        >>> motif_enrichment(foreground, background)
    {'UGAUUC': [0.72, 0.9302721574455114],
    'UAAACC': [2.4, 0.015461966703500418],
    'AAGUUACCU': [1.0, 0.0],
    'AAGCCUU': [1.1, 0.26809988836485743],
    'AGUUCUA': [0.44, 0.9053908005999413],
    'UUUCCCG': [1.65, 0.032318974952034396]}   
 """
    # make sure that the two dictionaries contain the same motifs
    for k in foreground:
        if k not in background.keys():
            background[k] = 0
    for k in background:
        if k not in foreground.keys():
            foreground[k] = 0

    # add pseudocount
    for k in foreground:
        foreground[k] += 1
    for k in background:
        background[k] += 1

    # calculate the total for each motif length
    total_foreground: Dict = dict()
    total_background: Dict = dict()
   
    for k in foreground:
        l = len(k)
        if l in total_foreground:
            total_foreground[l] += foreground[k]
        else:
             total_foreground[l] = foreground[k]
    
    for k in background:
        l = len(k)
        if l in total_background:
            total_background[l] += background[k]
        else:
             total_background[l] = background[k]


    # Main
    # Calculates enrichment
    enrichment_dict: Dict = dict()

    for k in foreground:

        enrichment_dict[k] = [((foreground[k] * total_background[len(k)]) /
                              (background[k] * total_foreground[len(k)]))]

    # Calculate p-value for foreground motifs
    # Create dictionary containing enrichment and p-values
    for k in foreground:
        n_bin = total_foreground[len(k)]
        p_bin = background[k] / total_background[len(k)]
        enrichment_dict[k].append( 1 - binom.cdf(foreground[k],
                                                       n_bin, p_bin))

    return enrichment_dict
