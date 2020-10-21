#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 10:10:35 2020

@author: Marit
"""

#conda install biopython 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs

# read sequence file
mydir= "/Users/Marit/Documents/programming/GitHub/htsinfer/tests/test_files/"
myfile= "SRR11971558_2.fastq"
myf=mydir+myfile
startlist=[]

for record in SeqIO.parse(myf, "fastq"):
    startlist.append(str(record.seq)) 

# create motif file
instances = [Seq("CC"), Seq('AA')]

m = motifs.create(instances)


for sequence in startlist:
    for pos, seq in m.instances.search(sequence):
       print("%i %s" % (pos, seq))



