##Issue description:

Given a set of core motifs (input: list of core motifs), check whether they represent
parts of longer motifs in sequenced reads (input: list of reads). 
Extend left and right if the nucleotide composition in these positions is very biased.
(output: list of extended motifs?)

That is, starting from all occurrences of a core motif in the reads, determine the 
frequencies of A,C,G,T at the left- and right-flanking positions, and calculate the
entropy of this distribution. Ideally, the motif can be extended if the flanking 
positions contain one specific nucleotide, i.e. entropy is 0. 

But sequencing errors could introduce some noise. So it may be necessary to introduce a cutoff.

Continue the extension greedily until the entropy becomes too high. 

##Proposal:

Input: list of core motifs & list of reads

Output: list of extended motifs

1) find core motif in every read (record start and end position)
2) determine the frequencies of nucleotides at the left- and right-flanking positions, 
compute the percentage of each nucleotide at the positions and the calculate 
the entropy at each position: -sum(P(i)*log2(P(i)) over all nucleotides(i)
3) if entropy<cutoff: extend motif in the right/left flanking position
4) continue until entropy>cutoff for both right and left flanking position 
(i.e. greedy approach)

