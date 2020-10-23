import time

test_sequences = ["CGTACTTGGTTCGATCAACCAGGCCGTAAAGTACGACGTCATCAAAAACGTTTAAACAAA",
                  "GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC",
                  "GCTGCTGGTCTTCATCCAACCTATGCTCGAACCATCGGTATTTCAGTTGATCATCGACGA",
                    "ACACGTCGTTATAACATGAAAGTACGTTCTGGACGCGGTTTTTCTTTGGATGAAATTCGT"]
#ATTUACGTAACGTGAAAAAGUAACTACGTCGACCTGTCAATGATGCTGATAUAACGAAATTTCAAGCTTTTCAAGCTUAA

test_not_list = 'GCTCAACGTGTTGCTCCTCGACCAGCCAAAGGTTCATTACGGCCAGTTGTTCGTGGTACC'

def count_motifs(all_sequences, min_motif_length, max_motif_length, nucleic_acid):
    """ Function that calculates the occurrence of all possible motifs in one or multiple sequences and
            returns a dictionary with all motifs within the specified length and their occurrence

        Args:
            input_sequences (list): list of sequences
            min_motif_length (int): minimal length of motif
            max_motif_length (int): maximal length of motif

        Returns:
            motif frequency (dict) with paired data {"motif_seq": frequency }

        Raises:
            TypeError: input_sequences is not a list
        """
    if not isinstance(all_sequences, list):
        print ('yay')
        raise TypeError('input_sequences is not a list')

    motif_freq = {}
    for i in range(0, len(all_sequences)):
        seq = all_sequences[i]
        user_choice = 0
        print (seq)
        if nucleic_acid.lower() == 'dna':
            valid_bases = 'ACTG'  # dict.fromkeys("012345789abcdef") might speed it up without hurting readability much.
        if nucleic_acid.lower() == 'rna':
            valid_bases = 'ACUG'
        if not all(base in valid_bases for base in seq):
            print('sequence', i+1,  'contains unknown characters')
            print('you have the following choices:')
            print('(1) ignore the non-nucleotide-characters and exclude them from the motifs')
            print('(2) include the non-nucleotide characters')
            print('(3) ignore the sequence with the non-nucleotide characters')
            user_choice = input('please specify using a number how you want to proceed: ')
        if user_choice == '3':
            print('sequence',i+1,'ignored')
            continue #skips this loop
        motifs_per_seq = []  # empty list for motives of this sequence
        for motif_length in range(min_motif_length, max_motif_length + 1):
            lst = []  # empty small list with kmers of length k
            for i in range(0, len(seq) - motif_length + 1):
                motif = seq[i: i + motif_length]
                #print (motif)
                if all(i in valid_bases for i in motif) or user_choice == '2':
                    if motif in motif_freq:
                        motif_freq[motif] += 1
                    else:
                        motif_freq[motif] = 1
                else:
                    print('not valid: ', motif)
                    # print (seq[i + motif_length])
            motifs_per_seq += lst  # add all kmers of length k to big list
    #return motif_freq    #motif_freq = count_motifs(motifs_total)  # create dictionary

    sorted_dict = dict(
        sorted(motif_freq.items(), key=lambda item: item[1], reverse=True))  # sorted dict in descending order
    #print(sorted_dict)
    return sorted_dict


if __name__ == "__main__":
    tic = time.perf_counter()
    #for i in range(0,100):
    count_motifs(all_sequences=test_sequences, min_motif_length=4, max_motif_length=5, nucleic_acid='dna')
    toc = time.perf_counter()
    print(f"Time to perform motif search: {toc - tic:0.4f} seconds")
