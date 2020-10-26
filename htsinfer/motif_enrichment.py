from timeit import timeit

foreground = {"TGATTC": 5, "TAAACC": 3, "AAGTTACCT": 1 "AAGCCTT": 1, "AGTTCTA": 0, "TTTCCCG": 5}
background = {"TGATTC": 3, "TAAACC": 5, "AAGCCTT": 0, "AGTTCTA": 1, "TTTCCCG": 5}

def motif_enrich(foreground_dic, background_dic):
	"""Calculates motif enrichment in two given dictionaries.

    Args:
        foreqorund_dic (dictionary) : .
        background_dic (dictionary) : .
        

    Returns:
        Dictionary including the motifs of the foreground dictionary and the enrichment 
        and p-value for each motif.

    Examples:
        
    """
    
    newDictForeground = dict()
    newDictEnrichForeground = dict()
    newDictBackground = dict()
    newDictEnrichBackground = dict()
    
    keysF = list(foregroundDic.keys()) #list containing all motifs
    keysB = list(backgroundDic.keys())
    lengthKeysF = [len(i) for i in foregroundDic.keys]  #list containing lengths of all motifs
    lengthKeysB = [len(i) for i in foregroundDic.keys]

    
    #foreground dictionary
    for i in keysF: 

    	if len(i) in lengthKeysF: 
	    	newDict[len(i)] += foregroundDic[i]
	    	#newDict[i] += foregroundDic[i] #creates dict with sum of counts for each motif with same lengt
	    else: 
	    	newDict[len(i)] = foregroundDic[i]
	    	#newDict[i] = foregroundDic[i]  #adds motif-lengths that are not yet there

	for i in newDict.values(): #calculates and adds enrichment score to newDict

		enrichScore = sqrt((len(foregroundDic) - i) / lengthF) 
		newDict[i] = enrichScore[i]

	for i in foregroundDic: 

		newDictEnrichForeground[i] = newDictForeground[len(i)].values()


	#background dictionary 
	for i in keysB: 

    	if len(i) in lengthKeysB: 
	    	newDictBackground[len(i)] += backgroundDic[i]
	    	#newDictBackground[i] += backgroundDic[i] #creates dict with sum of counts for each motif with same lengt
	    else: 
	    	newDict[len(i)] = backgroundDic[i]
	    	#newDictBackground[i] = backgroundDic[i]  #adds motif-lengths that are not yet there

	for i in newDict.values(): #calculates and adds enrichment score to newDict

		enrichScore = sqrt((len(backgroundDic) - i) / lengthF) 
		newDict[i] = enrichScore[i]

	for i in backgroundDic: 

		newDictEnrichBackground[i] = newDictBackground[len(i)].values() 

	#Calculating p-values 

	for i in foregroundDic: 
















	     

motif_enrich(foreground, background)

print("hello world")

