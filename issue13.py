test_sequence = "TTCTCCTAAGAATTGATTGATTCCCTTCATTTTCGGTAAGTACATGCTATATTATTTTATCCGCTTGAATCCTAGAACTTAGGACATCGATCGATCTCAATATTTGTTGATATAGGTAGTCAACAACAAAAGTCATCAAGATGGCTCCAAGCGGAAATAATATGATCCCAAATGGTCATTTTCACAAAGACTGGCAACGATACATTAAGACCTGGTTCAATCAACCAGCTCGTAAATACAGGCGTCACCAGAATCGTATTAAGAAAGCCAAAGAGTTGGCACCTAAACCAGCGGCGGGGCCTTTGCGACCAGTTGTACATTGCCCAACATTGCGATACCACACTAAAGTCCGTGCCGGACGTGGTTTTACTCTGCAAGAGCTGAAGGGTGCTGGTTTAAATAAACATTATGCCCGTACCA"

def find_motifs(seq, min_motif_length, max_motif_length):
    motif_freq = {}
    for motif_length in range(min_motif_length, max_motif_length + 1):
        print ("MOTIF LENGTH: ", motif_length)
        for i in range (0, len(seq) - motif_length + 1):
            kmer = seq[i : i + motif_length]
            print (kmer)

            if kmer in motif_freq:
                motif_freq[kmer] += 1
            else:
                motif_freq[kmer] = 1

    return motif_freq

if __name__ == "__main__":
    motif_dict = find_motifs(seq = test_sequence, min_motif_length = 3, max_motif_length = 6)
    print (motif_dict)