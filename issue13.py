
test_sequences = ["TTCCAGTCAGC", "CCGAGTCAGTAC", "TTCCAGTGGCCGA", "CCGGTAACTGGTAC"]

def find_motifs_per_seq(seq, min_motif_length, max_motif_length):
    motifs_per_seq = [] #empty big list for motives of this sequence
    for motif_length in range(min_motif_length, max_motif_length + 1):
        lst = [] # empty small list with kmers of length k
        for i in range (0, len(seq) - motif_length + 1):
            lst.append(seq[i : i + motif_length]) #kmer
        motifs_per_seq += lst #add all kmers of length k to big list
    return motifs_per_seq

def all_motifs(all_sequences, min_motif_length, max_motif_length):
    motifs_total = [] #empty big list for motives for all sequences
    for i in range(0,len(all_sequences)):
        motifs_per_seq = find_motifs_per_seq(all_sequences[i], min_motif_length, max_motif_length) #list with motives for one sequence
        motifs_total += motifs_per_seq
    return motifs_total

def count_motifs(motifs_total, min_motif_length, max_motif_length):
    motif_freq = {}
    # loop through the values in motifs_total and count them
    for motif in motifs_total:
        if motif in motif_freq:
            motif_freq[motif] += 1
        else:
            motif_freq[motif] = 1
    return motif_freq

def dict_motifs(all_sequences, min_motif_length, max_motif_length):
    motifs_total = all_motifs(all_sequences, min_motif_length, max_motif_length) # list with all motifs
    motif_freq = count_motifs(motifs_total, min_motif_length, max_motif_length) # create dictionary
    sorted_dict = dict(sorted(motif_freq.items(),key=lambda item: item[1],reverse=True)) # sorted dict in descending order
    print(sorted_dict)
    return sorted_dict

# test
dict_motifs(all_sequences = test_sequences, min_motif_length=4, max_motif_length=5)

#if __name__ == "__main__":
    #motif_dict = count_motifs(seq = test_sequence, min_motif_length = 3, max_motif_length = 6)
    #print (motif_dict)
            
