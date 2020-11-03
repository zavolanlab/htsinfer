#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 10:10:35 2020

@author: Marit
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join, getsize



#function
def motifHistogram(seq1, seq2, directory,nrbins):
    
    #-------------------------------------------------------------------------
    #create motif file of the 2 input motifs forward and revers
    instances = [Seq(seq1), Seq(seq2)]
    m = motifs.create(instances)
    r = m.reverse_complement()
    #-------------------------------------------------------------------------
    
    #get all reads of files in the folder
    files = [f for f in listdir(directory) if isfile(join(directory, f))] #list all files in the given directory

    #get only the fastqfiles in the folder
    onlyfastqfiles=[]
    for i in files:
        if i.endswith(".fastq"):
            onlyfastqfiles.append(i)
        else:
            print("other type of file discarded")
    
    #get only the filled fastq files in the list. This gives the list of fastq files ready to be checked for motifs
    filestoread=[]
    for i in onlyfastqfiles:
        if (getsize(directory+i)>0):
            filestoread.append(i)
        else:
            print("empty files are discareded")
            
    #-------------------------------------------------------------------------
    #loop over every fastq file
    #-------------------------------------------------------------------------
    for file in filestoread:
        #create an empty distance array where distances  between motifs will be stored for this fastqfile
        distance=[]
        
        #extract  the reads in the fastqfile, and store them in readlist
        readlist=[] #reset the readlist when a new file is opened
        for read in SeqIO.parse(directory+file, "fastq"):
            readlist.append(str(read.seq)) 
            
        #---------------------------------------------------------------------
        #loops over every sequence
        #---------------------------------------------------------------------
        for sequence in readlist: #takes the sequences from the input file
            motif1f=[]
            motif2f=[]
            motif1r=[]
            motif2r=[]
        
            # saves the positions in the forward motifs
            for pos, seq in m.instances.search(sequence): 
                if seq == m.instances[0]: #saves the positions of the first motif
                    motif1f.append(pos)
                if seq == m.instances[1]: #saves the positions of the second motif
                    motif2f.append(pos)
                    
            # saves the positions of the reverse motifs 
            for pos, seq in r.instances.search(sequence): 
                if seq == r.instances[0]: #saves the positions of the first motif
                    motif1r.append(pos)
                if seq == r.instances[1]: #saves the positions of the second motif
                    motif2r.append(pos)
     
            # calculate  distance between forward motif1 and motif2, and save it in distance array
            for i in motif1f:
                for j in motif2f:
                    distance.append(abs(i-j))
            
            
            # calculate distance between reverse motif1 and motif2, and save it in distance array
            for i in motif1r:
                for j in motif2r:
                    distance.append(abs(i-j))
                                        
    #---------------------------------------------------------------------
    plt.hist(distance, bins=10)
    plt.xlabel("distance (bp)")
    plt.ylabel("counts")
    plt.title("distance between given motifs")
    print('histogram with counts is created')
        


#-----------------------------------------------------------------------------
#import from frontend
#-----------------------------------------------------------------------------
mydir= "/Users/Marit/Documents/programming/GitHub/htsinfer/tests/test_files/"
nrbins=10
motifHistogram("CCG","GGG",mydir,nrbins)


 
     