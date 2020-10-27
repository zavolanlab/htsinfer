from timeit import timeit

foregroundDict = {"TGATTC": 5, "TAAACC": 3, "AAGTTACCT": 1, "AAGCCTT": 1, "AGTTCTA": 0, "TTTCCCG": 5}
backgroundDict = {"TGATTC": 3, "TAAACC": 5, "AAGCCTT": 0, "AGTTCTA": 1, "TTTCCCG": 5}

def motif_enrich(foreground, background):
	"""Calculates motif enrichment in two given dictionaries.

    Args:
        foreground (dictionary) : .
        background (dictionary) : .
        

    Returns:
        Dictionary including the motifs of the foreground dictionary and the enrichment 
        and p-value for each motif.

    Examples:
        
    """
	newDictForeground = dict()
	newDictEnrichForeground = dict()
	newDictBackground = dict()
	newDictEnrichBackground = dict()

    keysF = list(foreground.keys) #list containing all motifs
    keysB = list(background.keys())
    lengthKeysF = [len(i) for i in foreground.keys]  #list containing lengths of all motifs
    lengthKeysB = [len(i) for i in foreground.keys]

    
    #foreground dictionary
    for i in keysF: 

    	if len(i) in lengthKeysF: 
	    	newDict[len(i)] += foreground[i]
	    	#newDict[i] += foregroundDic[i] #creates dict with sum of counts for each motif with same lengt
	    else: 
	    	newDict[len(i)] = foreground[i]
	    	#newDict[i] = foregroundDic[i]  #adds motif-lengths that are not yet there

	for i in newDict.values(): #calculates and adds enrichment score to newDict

		enrichScore = sqrt((len(foreground) - i) / lengthF)
		newDict[i] = enrichScore[i]

	for i in foreground:

		newDictEnrichForeground[i] = newDictForeground[len(i)].values()


	#background dictionary 


	#Calculating p-values 

	for i in foregroundDic: 
















	     

motif_enrich(foregroundDict, backgroundDict)

print("hello world")

