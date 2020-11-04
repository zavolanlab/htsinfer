"""Calculate motif enrichment from sample data."""

import numpy as np
import scipy
from scipy.stats.distributions import binom


foreground = {"UGAUUC": 5, "UAAACC": 3, "AAGUUACCU": 1, "AAGCCUU": 1, "AGUUCUA": 1, "UUUCCCG": 5}
background = {"UGAUUC": 3, "UAAACC": 5, "AAGCCUUAU": 1, "AGUUCUA": 1, "UUUCCCG": 5, "UUGGAA": 7}


def motifEnrichment(
    foreground,
    background
): 
    """Calculates enrichment and p-values of motifs with similar lengths.

    Args:
        foreground: dictionary of motifs with counts
        background: dictionary of motifs with counts

    Returns:
        Dictionary including motifs, sum of counts of motifs with 
        similar lenght, enrichment score and p-value

    Example: 
        foreground = {"TGATTC": 5, "TAAACC": 3, "AAGTTACCT": 1, "AAGCCTT": 1, "AGTTCTA": 1, "TTTCCCG": 5}
        background = {"TGATTC": 3, "TAAACC": 5, "AAGCCTTAT": 1, "AGTTCTA": 1, "TTTCCCG": 5, "TTGGAA": 7}

        >>> motifEnrichment(foreground, background)
            {'TGATTC': [0.7333333333333334, 0.01646090534979423],
            'TAAACC': [0.7333333333333334, 0.21947873799725662],
            'AAGTTACCT': [1.375, 0.26337448559670795],
            'AAGCCTT': [1.6041666666666667, 0.26337448559670795],
            'AGTTCTA': [1.6041666666666667, 0.26337448559670795],
            'TTTCCCG': [1.6041666666666667, 0.01646090534979423]}
    """ 
   
     # Foreground dictionary
    newDictForeground = dict()
    newDictForegroundProb = dict()
    newDictForegroundFinal = dict()
    newDictForegroundBinom = dict()
    keyListF = list(foreground.keys())
    lengthListF = list()
    
    # Background dictionary
    newDictBackground = dict()
    newDictBackgroundProb = dict()
    keyListB = list(background.keys())
    lengthListB = list()
    
    # Foreground dictionary
    # Calulates sum of motifs of same length
    for i in keyListF: 
        
        if len(i) in lengthListF:
            newDictForeground[len(i)] += foreground[i]     
        else: 
            newDictForeground[len(i)] = foreground[i]
            lengthListF.append(len(i))
            
    # Calculates probability of occurrence of motifs
    for i in newDictForeground:
        
        probOfOccurrenceF = newDictForeground[i] / sum(foreground.values())
        newDictForegroundProb[i] = probOfOccurrenceF
    
    # Background dictionary
    # Calulates sum of motifs of same length
    for i in keyListB: 
        
        if len(i) in lengthListB:
            newDictBackground[len(i)] += background[i]     
        else: 
            newDictBackground[len(i)] = background[i]
            lengthListB.append(len(i))
    
    # Calculates probability of occurrence of motifs
    for i in newDictBackground:
        
        probOfOccurrenceF = newDictBackground[i] / sum(background.values())
        newDictBackgroundProb[i] = probOfOccurrenceF    
    
    # Main
    # Calculates enrichment 
    enrichmentDict = dict()
    
    for i in newDictForegroundProb:

        enrichmentDict[i] = newDictForegroundProb[i] / newDictBackgroundProb[i]    
    
    # Calculate p-value for foreground motifs    
    n = len(foreground)
    p = np.average(list(newDictBackgroundProb.values()))
        
    for i in foreground:
        
        r = foreground[i]
        newDictForegroundBinom[i] = binom.pmf(r, n, p)
       
    # Create dictionary containing enrichment and p-values 
    for i in foreground: 

        newDictForegroundFinal[i] = [enrichmentDict[len(i)], newDictForegroundBinom[i]]
    
    return(newDictForegroundFinal)


motifEnrichment(foreground, background)
