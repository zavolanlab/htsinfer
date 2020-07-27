##---------------------------------------------------------------
##-- The script extracts Gene ID's of a gene in an organism for all it's orthologues
##-- Example: Here we are extracting Homo Sapien's most abundant ribosomal 
##-- protein gene id's for all it's orthologues
##---------------------------------------------------------------

library(biomaRt)
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")

# List of Genes for which Gene ID's needs to be retrieved for it's orthologues 
genes = c("RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS5", 
            "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13")

## Prerequisite ****************
# Create an .csv file, before running the script, that will store the desired attributes needed for every organism 
df = read.csv('attributes.csv', header = FALSE)
df = sapply(df, as.character)
##******************************

for (gene in genes){
    
    # Creating empty dataframe that will store Gene ID's for all the orthologues for gene variable
    gene_ids_df = data.frame(matrix(nrow = 0, ncol=7))
    
    for (i in 1:nrow(df)){

        # Using getBM() to retrieve Gene ID for every organism in attributes.csv  
        id = getBM(filters = "external_gene_name",
                    attributes = c(df[i,"Type"], df[i,"Orthologue"], df[i,"Target %id"], 
                                    df[i,"Query %id"], df[i,"High Confidence"]),
                    values = gene,
                    mart = ensembl)
        id = unname(x)
        id = data.frame(df[i,"Species"], id, df[i,"Tax ID"])
        # print(x)
        # Merging dataframe for every organism's iteration
        gene_ids_df = rbind(as.matrix(gene_ids_df), as.matrix(id))
    }

    colnames(gene_ids_df) = c("Species", "Type", "Orthologue", "Target %id", 
                                "Query %id", "High Confidence", "Tax ID")
    file_name = paste("hsapiens_",gene, ".csv",sep="")
    write.csv(gene_ids_df, file=file_name, row.names = FALSE)
}
