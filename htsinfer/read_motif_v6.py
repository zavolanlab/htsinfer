# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 01:28:58 2020

@author: YB Moon (ymoon06)
"""
# assume motif <read
# use_N is for further implement as a parameter
def find_overlaps(motif, read, min_overlap, use_n=False):
    #check the type of arguments
    if not isinstance(motif, str) or\
        not isinstance(read, str) or\
        not isinstance(min_overlap, int) or\
        not isinstance(use_n, bool):
        raise TypeError('pass the right argument')
    # passing min_overlap <1 does not make sense
    if len(motif)==0 or min_overlap <1:
        raise ValueError('pass the right value')
        
    # compute partial overlaps of the motif at the start of the read (1st case)
    #list comprehension is faster. parallelizing is allowed
    partial_overlaps_start=[]
    if len(read)>= min_overlap:
        #doesnt make sense to do len(motif)+1 because then it would be the 2nd case
        #you start from min_overlap because over lap should satistfy certain threshold
        #As in motif[len(motif)-ov:], you dont need to specify the end, python understands this is until end 
        partial_overlaps_start= [ (0, ov/len(motif)) for ov in range(min_overlap, len(motif)) if read[0:ov]==motif[len(motif)-ov:]]
    
    # compute matches(overlaps) of the motif inside the read (2nd case)
    # you compare the full motif length
    full_overlaps=[]
    if len(read)>-len(motif):
        full_overlaps=[(i, 1) for i in range(0, len(read)-len(motif)+1) if read[i:i+len(motif)]==motif]
    
    # compute partial overlaps of the motif at the end of the read
    partial_overlaps_end=[]
    
    if len(read)>=min_overlap:
        partial_overlaps_end=[(len(read)-ov, ov/len(motif)) for ov in range(min_overlap, len(motif)) if read[len(read)-ov:]==motif[0:ov]]
    #print('success')
    # return the list of overlaps
    return partial_overlaps_start+ full_overlaps+ partial_overlaps_end
