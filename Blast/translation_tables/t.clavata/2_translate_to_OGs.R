library(rtracklayer)

# ----------------- connecting proteinID to assembly_IDs using blast results against the transcriptomic assembly ------

blast <- read.delim("filtered_blast_results.txt", stringsAsFactors = FALSE, header = FALSE)
table <- blast[, c(1,2)]
colnames(table) <- c("proteinID", "assemblysequenceID")

# ----------------- connecting to geneID using open access data annotation.txt (same source as scRNAseq data) -------

annotation <- read.delim("annotation.txt", stringsAsFactors = FALSE)
table$geneID <-NA

for (i in 1:length(table[,1])) {
  protein <- table[i,1]
  if (protein %in% annotation[,3]) {
    row <- which(annotation[,3] == protein)
    geneID <- annotation[row, 1]
    table[i,3] <- geneID
  }
}

# -------------unflattening matrix with OG (extracted clavata-column from filtered orthofinder.tsv) -----------
OGs <- read.table("orthogroups_IDs_filtered_ONLY_clavata.tsv", sep = "\t", header = TRUE, row.names = 1)
OGs_to_assemblysequences <- data.frame(assemblysequenceID = character(), orthogroup = character(), stringsAsFactors = FALSE)

for (og in rownames(OGs)) {
  if (OGs[og, 1] != "") {
    IDs <- unlist(strsplit(OGs[og, 1], ","))
    IDs <- trimws(IDs)
    # Create a temporary dataframe with the sub-elements and the rowname
    temp_df <- data.frame(assemblysequenceID = IDs, orthogroup = og, stringsAsFactors = FALSE)
    # Append the temporary dataframe to df_new
    OGs_to_assemblysequences <- rbind(OGs_to_assemblysequences, temp_df)
  }
}
# ------------- adding orthogroups -------------------------------------------
OGs <- OGs_to_assemblysequences
table$OG <- NA

for (i in 1:length(table[,1])) {
  assemblysequenceID <- table[i,2]
  if ((!is.na(assemblysequenceID)) && assemblysequenceID %in% OGs[,1]) {
    row <- which(OGs[,1] == assemblysequenceID)
    OG <- OGs[row, 2]
    table[i,4] <- OG
  }
}

write.table(table, file = "clavata_gene_to_protein_to_assemblysequenceID_to_OG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


