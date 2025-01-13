library(dplyr)
# ----------------- reading celltype annotation file (containing LASC_IDs) and cleaning it  ------

table <- read.csv("Markers.SCT.integrated.majorCluster.all.2024.csv", header = TRUE, row.names = NULL)
table <- table[, c("p_val_adj", "cluster", "lascID")]
table <- table[order(table$cluster, table$p_val_adj), ] # making sure that it's sorted the right way
table <- as.data.frame(table %>% 
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 50) %>% 	# subsetting top 50 genes per cluster
  ungroup())

# -------------unflattening matrix with OG (extracted lasc-column from filtered orthofinder.tsv) -----------
OGs <- read.table("lasc_OGs.tsv", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
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

write.table(table, file = "LASC_id_to_OG.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


