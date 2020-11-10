#!/usr/bin/env python
#This function searches a sequence file with a set of motifs to produce a list of positions where the motif is present in each sequence
    
import re
from Bio import SeqIO   
import sys

sys.stdout=open("output_positions.txt","w")


def reader ():
    motif = "ATAC" #This should also be read from a file
    for record in SeqIO.parse("transcripts.fasta", "fasta"):
        print (record.id)
        sequence = str(record.seq)
        s = re.finditer(motif,sequence)

        for match in s:
            print (match.start())

sys.stdout.close()

import sys
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sys.stdout=open("small_hist.txt","w")

#This function draws a histogram of the position frequencies for each motif
    
with open('output_positions.txt', 'r') as fd:
    lines = fd.read().split()
    counter = Counter(lines)
    # sorts items
    items = sorted(counter.items(), key=lambda x: int(x[0]))
    # prints desired output
    for k, repetitions in items:
        print (k,'\t', repetitions)
    ind = np.arange(len(repetitions))
    plt.bar(k, ind, width=0.6, color='b')
    plt.show()
 
sys.stdout.close()