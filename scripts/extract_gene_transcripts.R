##---------------------------------------------------------------
##-- The script extracts transcripts of a gene in an organism for all it's orthologues
##-- Example: Here we are extracting Saccharomyces cerevisiae's most abundant ribsomal 
##-- protein gene transcripts for all it's orthologues
##---------------------------------------------------------------

library(biomaRt)
library(seqinr)

# List of Genes for transcripts needs to be retrieved
genes_fungi = c("RPS0A", "RPS0B", "RPS2", "RPS3", "RPS1A",
                    "RPS1B", "RPS4A", "RPS4B", "RPS5", "RPS6A",
                        "RPS6B", "RPS7A", "RPS7B", "RPS8A", "RPS8B")
number_genes = length(genes_fungi)

# Create .fasta file for each respective gene 
for (i in 1:number_genes){
    gene_transcript_file = paste("scerevisiae_", genes_fungi[i], ".fasta", sep = "")
    cat( "", file = gene_transcript_file, append = FALSE, sep = "")
}

# Looping for every Gene in genes_fungi
for (i in 1:number_genes){
    
    # Reading .csv file that was created by extract_gene_id script
    gene_id_file = paste("scerevisiae_", genes_fungi[i], ".csv", sep = "" )
    df = read.csv(gene_id_file, header = TRUE )
    df = sapply(df, as.character)
    
    # Iterating for every organism in .csv file
    for(j in 1:nrow(df)){

        organism_name = df[j,"Species"]
        gene_ensembl_name = df[j,"Orthologue"]

        # If no Gene ID for that organism then skip it 
        if(is.na(gene_ensembl_name) || gene_ensembl_name="")
        {
            next
        }

        # Here dataset refers to organism's dataset in ensembl mart, example: aclavatus_eg_gene
        dataset = paste(organism_name, "_eg_gene", sep = "")

        # For other Ensembls like Fungi, Plants, Bacteria
        ensembl = useMart(host = "https://fungi.ensembl.org", 
                            biomart = "fungi_mart", dataset = dataset, 
                            port = 443)
        
        # For regular Ensembl Biomart 
        # ensembl <- useEnsembl(biomart = "ensembl", 
        #                        dataset = dataset, 
        #                        mirror = "asia")
        
        # Dataframe with Transcript ID and cdna sequence
        transcript_df = biomaRt::getBM(filters = "ensembl_gene_id",
                           attributes = c("ensembl_transcript_id", "cdna"),  
                           values = gene_ensembl_name,
                           mart = ensembl)
        
        gene_transcript_file = paste("scerevisiae_", genes_fungi[i], ".fasta", sep = "")
        for (k in 1:nrow(transcript_df)){
            
            # print(transcript_df[k,1])
        
            id = paste(genes_fungi[i], gene_ensembl_name, 
                        transcript_df[k,1], organism_name, 
                            df[j,"Tax ID"], sep="|")
            # Appending transcripts to the gene's fasta file
            write.fasta(sequences = transcript_df[k,2], names=id, 
                                file.out = gene_transcript_file, 
                                open = "a", nbchar = 60, 
                                as.string = TRUE)
        
            }
        
    }
    
}

## Post-run ****************
## Once this is done for all the genes, we need to run getBM() one for time 
## to get the transcripts for the organism itself, as in the case we are running
## the script for Scerevisiae so we will run the getBM() for scerevisiae to 
## retrieve transcripts for all the aformentioned genes.
##******************************