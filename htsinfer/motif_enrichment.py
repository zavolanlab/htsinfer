"""Calculate motif enrichment from sample data."""

from timeit import timeit
import scipy
from scipy.stats.distributions import binom
import sys


foreground = {"TGATTC": 5, "TAAACC": 3, "AAGTTACCT": 1, "AAGCCTT": 1, "AGTTCTA": 0, "TTTCCCG": 5}
background = {"TGATTC": 3, "TAAACC": 5, "AAGCCTTAT": 0, "AGTTCTA": 1, "TTTCCCG": 5, "TTGGAA": 7}


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

    Raises:
        .
    """ 
   
     # Foreground dictionary
    newDictForeground = dict()
    newDictForegroundProb = dict()
    newDictFinalF = dict()
    keyListF = list(foreground.keys())
    lengthListF = list()
    
    # Background dictionary
    newDictBackground = dict()
    newDictBackgroundProb = dict()
    newDictFinalB= dict()
    keyListB = list(background.keys())
    lengthListB = list()
    
    # Foreground dictionary
    # Calulates sum of motifs of same length
    for i in keyListF: 
        
        if len(i) in lengthListF:
            newDictForeground[len(i)] += 1      
        else: 
            newDictForeground[len(i)] = 1
            lengthListF.append(len(i))
            
    # Calculates probability of occurrence of motifs
    for i in newDictForeground:
        
        probOfOccurrenceF = newDictForeground[i] / len(foreground)
        newDictForegroundProb[i] = probOfOccurrenceF
      
    # Background dictionary
    # Calulates sum of motifs of same length
    for i in keyListB: 
        
        if len(i) in lengthListB:
            newDictBackground[len(i)] += 1      
        else: 
            newDictBackground[len(i)] = 1
            lengthListB.append(len(i))
            
    # Calculates probability of occurrence of motifs
    for i in newDictBackground:
        
        probOfOccurrenceF = newDictBackground[i] / len(background)
        newDictBackgroundProb[i] = probOfOccurrenceF
    
        
    # Calculates enrichment 
    enrichmentDict = dict()
    
    for i in newDictForegroundProb:

        enrichmentDict[i] = newDictForegroundProb[i] / newDictBackgroundProb[i]
    
    # Adds enrichment values
        
    for i in foreground: 

        newDictFinalF[i] = enrichmentDict[len(i)]
    
    for i in background: 

        newDictFinalB[i] = enrichmentDict[len(i)]
        
    print(newDictFinalF)
    print(newDictFinalB)

motifEnrichment(foreground, background)

