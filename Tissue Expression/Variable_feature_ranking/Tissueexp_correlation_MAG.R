library(ggplot2)

LASC_data <- read.delim("C:/Users/Cleiti/Desktop/Seurat_spiders/bulkRNA_pertissue/summed_lasc.expression.txt", header = TRUE, sep = "\t", row.names = 1)
clavata_data <- read.delim("C:/Users/Cleiti/Desktop/Seurat_spiders/bulkRNA_pertissue/summed_clavata_expression.txt", header = TRUE, sep = "\t", row.names = 1)

LASC_data_MAG <- LASC_data[, c("Duct", "Sac", "Tail")]
clavata_data_MAG <- clavata_data[, c("Duct", "Sac", "Tail")]
# filtering out rows with too low tpm values (if tpm in any column is <= 10) = assumed noise
rows_to_remove <- c()
for (i in 1:nrow(LASC_data_MAG)) { 
  if (!(all(LASC_data_MAG[i, ] > 10))) {
    rows_to_remove <- c(rows_to_remove, rownames(LASC_data_MAG)[i])
  }
}
LASC_data_MAG_filtered <- LASC_data_MAG[!rownames(LASC_data_MAG) %in% rows_to_remove, ]
rows_to_remove <- c()
for (i in 1:nrow(clavata_data_MAG)) { 
  if (!(all(clavata_data_MAG[i, ] > 10))) {
    rows_to_remove <- c(rows_to_remove, rownames(clavata_data_MAG)[i])
  }
}
clavata_data_MAG_filtered <- clavata_data_MAG[!rownames(clavata_data_MAG) %in% rows_to_remove, ]
# log normalization (natural log)
log_LASC_data_MAG <- log(LASC_data_MAG_filtered)
log_clavata_data_MAG <- log(clavata_data_MAG_filtered)
# filtering out orthogroups not present in both datasets
common_rows <- intersect(rownames(log_LASC_data_MAG), rownames(log_clavata_data_MAG))
log_LASC_data_MAG_filtered <- log_LASC_data_MAG[common_rows, , drop = FALSE]
log_clavata_data_MAG_filtered <- log_clavata_data_MAG[common_rows, , drop = FALSE]
# calculating feature variances
feature_variances_lasc <- apply(log_LASC_data_MAG_filtered, 1, var, na.rm=TRUE)
feature_variances_clavata <- apply(log_clavata_data_MAG_filtered, 1, var, na.rm=TRUE)
sorted_variances_lasc <- sort(feature_variances_lasc, decreasing = TRUE)
sorted_variances_clavata <- sort(feature_variances_clavata, decreasing = TRUE)
sorted_variances_df_lasc <- data.frame(orthogroup = names(sorted_variances_lasc), rank_LASC = seq_along(sorted_variances_lasc))
sorted_variances_df_clavata <- data.frame(orthogroup = names(sorted_variances_clavata), rank_clavata = seq_along(sorted_variances_clavata))

# Correlation plot
merged_df <- merge(sorted_variances_df_lasc, sorted_variances_df_clavata, by = "orthogroup")
correlation_plot <- ggplot(merged_df, aes(x = rank_LASC, y = rank_clavata)) +
  geom_point() +                 
  geom_smooth(method = "lm") +   
  labs(x = "rank_lasc", y = "rank_clavata", 
       title = "Correlation of Orthogroup Variance Ranking in Major Ampullate Gland sections") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))
print(correlation_plot)
model <- lm(rank_clavata ~ rank_LASC, data = merged_df)
slope <- coef(model)[2]
print(paste("Correlation value: ", slope))
