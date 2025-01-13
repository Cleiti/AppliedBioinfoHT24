library(ggplot2)

LASC_data <- read.delim("C:/Users/Cleiti/Desktop/Seurat_spiders/bulkRNA_pertissue/summed_lasc.expression.txt", header = TRUE, sep = "\t", row.names = 1)
clavata_data <- read.delim("C:/Users/Cleiti/Desktop/Seurat_spiders/bulkRNA_pertissue/summed_clavata_expression.txt", header = TRUE, sep = "\t", row.names = 1)

# ------------------ HEATMAP: OG VS TISSUES (top OGs by variance)

# BRIDGE glands data prep
LASC_data_glands <- LASC_data[, c("Aggregate.gland", "Major.ampullate.gland")]
rows_to_remove <- c()
for (i in 1:nrow(LASC_data_glands)) { # filtering out rows with too low tpm values (most likely noise)
  if (!(all(LASC_data_glands[i, ] > 10))) {
    rows_to_remove <- c(rows_to_remove, rownames(LASC_data_glands)[i])
  }
}
LASC_data_glands_filtered <- LASC_data_glands[!rownames(LASC_data_glands) %in% rows_to_remove, ]
log_LASC_data_glands <- log(LASC_data_glands_filtered)
feature_variances <- apply(log_LASC_data_glands, 1, var, na.rm=TRUE)
sorted_vars <- sort(feature_variances, decreasing = TRUE)
selected <- sorted_vars[1:30]
selected_IDs <- names(selected)
selected_rownumbers <- which(rownames(log_LASC_data_glands) %in% selected_IDs)
log_LASC_data_glands_selected <- log_LASC_data_glands[selected_rownumbers, , drop = FALSE]
log_LASC_data_glands_selected$Row <- rownames(log_LASC_data_glands_selected)
log_LASC_data_glands_long <- reshape2::melt(log_LASC_data_glands_selected, id.vars = "Row", variable.name = "Tissue", value.name = "Value")

# BRIDGE only Major Ampullate Gland data prep
LASC_data_MAG <- LASC_data[, c("Duct", "Sac", "Tail")]
rows_to_remove <- c()
for (i in 1:nrow(LASC_data_MAG)) { # filtering out rows with too low tpm values (most likely noise)
  if (!(all(LASC_data_MAG[i, ] > 10))) {
    rows_to_remove <- c(rows_to_remove, rownames(LASC_data_MAG)[i])
  }
}
LASC_data_MAG_filtered <- LASC_data_MAG[!rownames(LASC_data_MAG) %in% rows_to_remove, ]
log_LASC_data_MAG <- log(LASC_data_MAG_filtered)
feature_variances <- apply(log_LASC_data_MAG, 1, var, na.rm=TRUE)
sorted_vars <- sort(feature_variances, decreasing = TRUE)
selected <- sorted_vars[1:30]
selected_IDs <- names(selected)
selected_rownumbers <- which(rownames(log_LASC_data_MAG) %in% selected_IDs)
log_LASC_data_MAG_selected <- log_LASC_data_MAG[selected_rownumbers, , drop = FALSE]
log_LASC_data_MAG_selected$Row <- rownames(log_LASC_data_MAG_selected)
log_LASC_data_MAG_long <- reshape2::melt(log_LASC_data_MAG_selected, id.vars = "Row", variable.name = "Tissue", value.name = "Value")

# histogram over log(tpm) values (after tpm filtering)
library(tidyr)
log_LASC_data_glands_forhistogram <- pivot_longer(log_LASC_data_glands, cols = everything(), names_to = "Column", values_to = "Value")
ggplot(log_LASC_data_glands_forhistogram, aes(x = Value, fill = Column)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, color = "black", position = "identity") +
  facet_wrap(~Column, scales = "free_y") +
  labs(title = "Value Distributions by Column", x = "Values", y = "Frequency") +
  theme_minimal()
log_LASC_data_MAG_forhistogram <- pivot_longer(log_LASC_data_MAG, cols = everything(), names_to = "Column", values_to = "Value")
ggplot(log_LASC_data_MAG_forhistogram, aes(x = Value, fill = Column)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, color = "black", position = "identity") +
  facet_wrap(~Column, scales = "free_y") +
  labs(title = "Value Distributions by Column", x = "Values", y = "Frequency") +
  theme_minimal()

# heatmap BRIDGE glands
ggplot(log_LASC_data_glands_long, aes(Tissue, Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Tissue", y = "Orthogroup", fill = "Value") +
  ggtitle("L.Sclopetarius most variable OGs within glands (lognormalized)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))

# heatmap BRIDGE only Major Ampullate Gland
ggplot(log_LASC_data_MAG_long, aes(Tissue, Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Tissue", y = "Orthogroup", fill = "Value") +
  ggtitle("L.Sclopetarius most variable OGs within major ampullate gland sections (lognormalized)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))

# CLAVATA glands data prep
clavata_data_glands <- clavata_data[, c("Aggregate.gland", "Major.ampullate.gland")]
rows_to_remove <- c()
for (i in 1:nrow(clavata_data_glands)) { # filtering out rows with too low tpm values
  if (!(all(clavata_data_glands[i, ] > 10))) {
    rows_to_remove <- c(rows_to_remove, rownames(clavata_data_glands)[i])
  }
}
clavata_data_glands_filtered <- clavata_data_glands[!rownames(clavata_data_glands) %in% rows_to_remove, ]
log_clavata_data_glands <- log(clavata_data_glands_filtered)
feature_variances <- apply(log_clavata_data_glands, 1, var, na.rm=TRUE)
sorted_vars <- sort(feature_variances, decreasing = TRUE)
selected <- sorted_vars[1:30]
selected_IDs <- names(selected)
selected_rownumbers <- which(rownames(log_clavata_data_glands) %in% selected_IDs)
log_clavata_data_glands_selected <- log_clavata_data_glands[selected_rownumbers, , drop = FALSE]
log_clavata_data_glands_selected$Row <- rownames(log_clavata_data_glands_selected)
log_clavata_data_glands_long <- reshape2::melt(log_clavata_data_glands_selected, id.vars = "Row", variable.name = "Tissue", value.name = "Value")

# CLAVATA only Major Ampullate Gland data prep
clavata_data_MAG <- clavata_data[, c("Duct", "Sac", "Tail")]
rows_to_remove <- c()
for (i in 1:nrow(clavata_data_MAG)) { # filtering out rows with too low tpm values
  if (!(all(clavata_data_MAG[i, ] > 10))) {
    rows_to_remove <- c(rows_to_remove, rownames(clavata_data_MAG)[i])
  }
}
clavata_data_MAG_filtered <- clavata_data_MAG[!rownames(clavata_data_MAG) %in% rows_to_remove, ]
log_clavata_data_MAG <- log(clavata_data_MAG_filtered)
feature_variances <- apply(log_clavata_data_MAG, 1, var, na.rm=TRUE)
sorted_vars <- sort(feature_variances, decreasing = TRUE)
selected <- sorted_vars[1:30]
selected_IDs <- names(selected)
selected_rownumbers <- which(rownames(log_clavata_data_MAG) %in% selected_IDs)
log_clavata_data_MAG_selected <- log_clavata_data_MAG[selected_rownumbers, , drop = FALSE]
log_clavata_data_MAG_selected$Row <- rownames(log_clavata_data_MAG_selected)
log_clavata_data_MAG_long <- reshape2::melt(log_clavata_data_MAG_selected, id.vars = "Row", variable.name = "Tissue", value.name = "Value")

# bridge only MAG gland just compare OG50 values:
print(log_clavata_data_MAG["OG0000050",])

# histogram over log(tpm) values (after tpm filtering)
library(tidyr)
log_clavata_data_glands_forhistogram <- pivot_longer(log_clavata_data_glands, cols = everything(), names_to = "Column", values_to = "Value")
ggplot(log_clavata_data_glands_forhistogram, aes(x = Value, fill = Column)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, color = "black", position = "identity") +
  facet_wrap(~Column, scales = "free_y") +
  labs(title = "Value Distributions by Column", x = "Values", y = "Frequency") +
  theme_minimal()
log_clavata_data_MAG_forhistogram <- pivot_longer(log_clavata_data_MAG, cols = everything(), names_to = "Column", values_to = "Value")
ggplot(log_clavata_data_MAG_forhistogram, aes(x = Value, fill = Column)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, color = "black", position = "identity") +
  facet_wrap(~Column, scales = "free_y") +
  labs(title = "Value Distributions by Column", x = "Values", y = "Frequency") +
  theme_minimal()

# heatmap CLAVATA glands
ggplot(log_clavata_data_glands_long, aes(Tissue, Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Tissue", y = "Orthogroup", fill = "Value") +
  ggtitle("T.Clavata most variable OGs within glands (lognormalized)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))

# heatmap CLAVATA only Major Ampullate Gland
ggplot(log_clavata_data_MAG_long, aes(Tissue, Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Tissue", y = "Orthogroup", fill = "Value") +
  ggtitle("T.Clavata most variable OGs within major ampullate gland sections (lognormalized)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))
