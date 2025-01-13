library(rtracklayer)

# ------------------ connecting gene_id to protein_id using NepCla refgenome .gtf file --------------
gtf_data <- import("genomic.gtf", format = "gtf")
gtf_df <- as.data.frame(gtf_data)

subset_gtf <- gtf_df[, c("gene_id", "protein_id")]
unique_only <- dplyr::distinct(subset_gtf)
cleaned_df <- unique_only[complete.cases(unique_only[, c("gene_id", "protein_id")]), ]

# adjusting gene names to converge with scRNAseq dataset
cleaned_df$gene_id <- gsub("_", "-", cleaned_df$gene_id)
cleaned_df$gene_id <- paste0(cleaned_df$gene_id, "-RA")

# ----------------- connecting to assembly_IDs using blast results against the transcriptomic assembly ------

table <- cleaned_df 
annotation <- read.delim("filtered_blast_results.txt", stringsAsFactors = FALSE, header = FALSE)
annotation$V1 <- sapply(annotation$V1, function(x) strsplit(x, "_")[[1]][3]) # fixing format to align with proteinIDs

table$assemblysequenceID <-NA

for (i in 1:length(table[,1])) {
  protein <- table[i,2]
  if (protein %in% annotation[,1]) {
    row <- which(annotation[,1] == protein)
    assemblysequenceID <- annotation[row, 2]
    table[i,3] <- assemblysequenceID
  }
}

# -------------unflattening matrix with OG (extracted clavipes-column from filtered orthofinder.tsv) -----------
OGs <- read.table("OGs_only_clavipes.tsv", sep = "\t", header = TRUE)
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
  assemblysequenceID <- table[i,3]
  if ((!is.na(assemblysequenceID)) && assemblysequenceID %in% OGs[,1]) {
    row <- which(OGs[,1] == assemblysequenceID)
    OG <- OGs[row, 2]
    table[i,4] <- OG
  }
}

write.table(table, file = "clavipes_gene_to_protein_to_assemblysequenceID_to_OG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


