from Bio import SeqIO
import numpy as np
from tqdm import tqdm
import pandas as pd
import math


def count_overlapping_motif_in_read(motif, read, overlap_cuttoff_ratio ):
    """ Given a motif and a read, compare the motif against all substrings of the read to search for (partial) overlaps.
    
    Args:
        motif: a short string, e.g., 'TTCTCCTAAGAATTGATTGAT'
        read: a long string, whose length needs to be greater than that of the motif, e.g., 
                    'TTCTCCTTTCTCCTAAGAATTGATTGATTCTCCTAAGAATTGATTGATAAGAATTGATTGAT'
        overlap_cuttoff_ratio: a float between 0 and 1. Only overlaps larger than len(motif) * overlap_cuttoff_ratio will be taken into account.
    Returns:
        motif_overlap_count_dict: a dictionary that receives the number of partially aligned characters as input produce the number of substrings 
             in the read that has possess this amount of partially aligned characters as output when compared to the given motif; 
             
        overlaps_dict: collects all partially aligned substrings for a given number of partially aligned characters.
    """
    len_motif = len(motif)


    if len_motif > len(read):
        print('read: ', read)
        print('motif: ', motif)
        raise ValueError('the motif length needs to be shorter than the read length!') 
    
    cutoff_length = math.ceil(len_motif * overlap_cuttoff_ratio)
    motif_overlap_count_dict = {x: 0 for x in range(cutoff_length, len_motif + 1)}
    overlaps_dict = {x: [] for x in range(cutoff_length, len_motif + 1)}
    
    for i in range(len(read) - len_motif + 1):

        # extract a sub-string of the read of the length of the motif. 
        sub_read = read[i:i+ len_motif ]

        # return an array of True-False value -- matched character corresponds to True, else False.
        motif_matching_bool = ( np.array(list(motif) ) == np.array(list(sub_read))  )

        # count the number of matched characters.
        motif_matching_count = np.sum( motif_matching_bool.astype(int) )

        if motif_matching_count >= cutoff_length:
            # update the count.
            motif_overlap_count_dict[motif_matching_count] += 1

            # store the sub-string of the read.
            overlaps_dict[motif_matching_count].append(sub_read)

    return motif_overlap_count_dict, overlaps_dict
    
def analyze_single_motif_occurence(motif, overlap_cuttoff_ratio, fasta_path = './transcripts.fasta'):
    """ Given a motif and a fasta file, compare the motif against all reads in the file.
    Args:
        motif: a short string, e.g., 'TTCTCCTAAGAATTGATTGAT'
        overlap_cuttoff_ratio: a float between 0 and 1. Only overlaps larger than len(motif) * overlap_cuttoff_ratio will be taken into account.
        fasta_path: the path to the fasta file.
    Returns:
        motif_count_over_reads_dict and overlaps_across_reads_dict assemble outputs of the function count_overlapping_motif_in_read for all reads. 
    """
    all_reads = list( SeqIO.parse(fasta_path, 'fasta') )
    all_reads = [read for read in all_reads if str(read.seq ) != 'Sequenceunavailable' ]


    motif_count_over_reads_dict = {}
    overlaps_across_reads_dict = {}

    for read_count in tqdm(range(len(all_reads))):
        read = str( all_reads[read_count].seq )
        read_id = all_reads[read_count].id
        motif_overlap_count_dict, overlaps_dict = count_overlapping_motif_in_read(motif, read, overlap_cuttoff_ratio)
        motif_count_over_reads_dict[read_id] = motif_overlap_count_dict
        overlaps_across_reads_dict[read_id] = overlaps_dict
        
    return motif_count_over_reads_dict, overlaps_across_reads_dict

