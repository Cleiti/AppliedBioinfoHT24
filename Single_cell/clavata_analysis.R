library(Seurat) #v5
library(DoubletFinder)
library(clustree)
library(patchwork)
library(ggplot2)
library(dplyr)


data <- Read10X(data.dir = "C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_data") # folder with 3 files
obj <- CreateSeuratObject(counts = data, project = "scRNA_t.clavata")

plot_QC <- function(seurat_obj) {
  
  # violin plots
  print(((VlnPlot(seurat_obj,  features = c("nFeature_RNA"), pt.size=0, y.max=4000, raster=FALSE) + ylab("Genes per cell") + xlab(NULL) + ggtitle(NULL) +theme(legend.position = "none")) +
      (VlnPlot(seurat_obj,  features = c("nCount_RNA"), pt.size=0, y.max=5000, raster=FALSE) + ylab("UMI per cell") + xlab(NULL)+ ggtitle(NULL) +theme(legend.position = "none")) +
      plot_layout(ncol = 2) +  theme(axis.title = element_text(size=10), axis.text = element_text(size=8)) + NoLegend()))
  
  # same info but in distribution plots
  print(seurat_obj@meta.data %>% 
    ggplot(aes(x=nCount_RNA)) +  # UMI
    geom_density(alpha = 0.2) + 
    theme_classic() +
    ylab("log10 cell density") +
    scale_x_log10())
  
  print(seurat_obj@meta.data %>% 
    ggplot(aes(x=nFeature_RNA)) + #Genes
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10())
}

#initial raw data QC:
plot_QC(obj)
  # (compare with normal: https://github.com/hbctraining/scRNA-seq_online/blob/2b6165fd3867a9cbe242bbf5d1699c4b05ba70cf/lessons/04_SC_quality_control.md)

#---------------------------------------------------------------------
# filtering by UMI:
obj_filtered <- subset(obj, subset = nCount_RNA <= 10000)
plot_QC(obj_filtered)

# filtering by genes (to avoid taking in second maximum on distribution plot - putative doublets)
obj_filtered <- subset(obj, subset = nFeature_RNA < 1500)
plot_QC(obj_filtered)

#---------------------------------------------------------------------
# filtering by genes only present in filtered orthofinder
annotation <- read.table("C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_gene_to_protein_to_assemblysequenceID_to_OG.tsv", sep = "\t", header = TRUE)

features_filtered_out <- c()
for (geneID in rownames(obj_filtered)) {
  if (!(geneID %in% annotation[,3])) {
    features_filtered_out <- cbind(features_filtered_out, geneID)
  }
}
obj_filtered_by_OG <- subset(obj_filtered, features = setdiff(rownames(obj_filtered), features_filtered_out))
# -> 1/4 of features got removed

plot_QC(obj_filtered_by_OG)

# ANALYSIS CONTINUED: ---------------------------------------------DoubletFinder

# OBS: it should run PER SAMPLE (in a loop). 
# The source paper says they had 3 spiders = 3 samples
# BUT the data has no sample information.
# -> running on the whole thing

# Running SCTransform, PCA
obj_SCT <- SCTransform(obj_filtered_by_OG, vars.to.regress= "nFeature_RNA")  #  normalization
obj_PCA <- RunPCA(obj_SCT) # dimensionality reduction

# Calculating significant PCs
#(source: https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html)
stdv <- obj_PCA[["pca"]]@stdev
sum.stdv <- sum(obj_PCA[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2) 

# also visualized
ElbowPlot(obj_PCA, ndims = 40) + geom_vline(xintercept = 12, linetype="dashed", color = "green") 

# Running UMAP, FindNeighbors, FindClusters
obj <- RunUMAP(obj_PCA, reduction = "pca", dims = 1:pcs)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:pcs)
obj <- FindClusters(obj, resolution = 0.1)

# the next doubletfinder steps are identical to tutorials, we don't have to understand every single one

# Running PK identification (no ground-truth)
sweep.list <- paramSweep(obj, PCs = 1:pcs, num.cores = 1, sct=T) # more than 1 is pain in the ass on windows
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats) 

# Calculating optimal pk
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

# Estimating homotypic doublet proportion
annotations <- obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(optimal.pk * nrow(obj@meta.data))
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# Running DoubletFinder
obj <- doubletFinder(seu = obj, PCs = 1:pcs, pK = optimal.pk, nExp = nExp.poi.adj, sct=T)

# Subsetting singlets
df_classifications_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
obj_singlets <- subset(obj, subset = !!sym(df_classifications_col) == "Singlet")
# -> only filtered out ~50 cells
obj <- obj_singlets

# filtering results
plot_QC(obj)
#-------------------------------------------------------------------------downstream analysis

# Rerunning SCTransform, PCA (less data after doubletfinder -> redoing these even though we did them before)
obj <- SCTransform(obj, vars.to.regress= "nFeature_RNA")  
obj <- RunPCA(obj)

# Calculating significant PCs
#(source: https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html)
stdv <- obj[["pca"]]@stdev
sum.stdv <- sum(obj[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

# also visualized
ElbowPlot(obj, ndims = 40) + geom_vline(xintercept = 12, linetype="dashed", color = "green") 

# Running UMAP, FindNeighbors, FindClusters
obj <- RunUMAP(obj, reduction = "pca", dims = 1:pcs)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:pcs)
obj <- FindClusters(obj, resolution = seq(0.1, 0.4, by=0.1))

# how stable are the clusters calculated at different resolutions?
# (= how sure is the software that these celltypes are DISTINCT?)
clustree(obj@meta.data[,grep("SCT_snn_res.", colnames(obj@meta.data))],prefix = "SCT_snn_res.")

(DimPlot(obj,
         group.by=grep("SCT_snn_res", colnames(obj@meta.data), value = TRUE)[1:4],
         shuffle = TRUE,
         ncol=2 ,
         reduction = "umap",
         label = T) + plot_annotation(title = 'Clustering resolutions')) & theme(legend.key.size = unit(0.2, 'cm'))

saveRDS(obj, "C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_results/clavata_analyzed.rds")
#-----------------------findmarkers = significant features for every cluster

# running findallmarkers on all resolutions
DefaultAssay(object = obj) <- "SCT"

Idents(object = obj) <- "SCT_snn_res.0.1"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.1 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.1")
Idents(object = obj) <- "SCT_snn_res.0.2"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.2 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.2")
Idents(object = obj) <- "SCT_snn_res.0.3"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.3 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.3")
Idents(object = obj) <- "SCT_snn_res.0.4"
obj <- PrepSCTFindMarkers(obj)
markers_res.0.4 <- FindAllMarkers(obj, group.by = "SCT_snn_res.0.4")

write.csv(markers_res.0.1, file = "C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_results/markers_res.0.1.csv", row.names = TRUE)
write.csv(markers_res.0.2, file = "C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_results/markers_res.0.2.csv", row.names = TRUE)
write.csv(markers_res.0.3, file = "C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_results/markers_res.0.3.csv", row.names = TRUE)
write.csv(markers_res.0.4, file = "C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_results/markers_res.0.4.csv", row.names = TRUE)


#--------------------- annotation via OG comparison to LASC celltypes
annotation <- read.table("C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_gene_to_protein_to_assemblysequenceID_to_OG.tsv", sep = "\t", header = TRUE)
OGs <- annotation # for clavata
OGs <- OGs[, c(3,4)]

# prepare lasc cluster information
lasc_annotation <- read.table("C:/Users/Cleiti/Desktop/Seurat_spiders/translation_tables/LASC_id_to_OG.tsv", sep = "\t", header = TRUE)
lasc_OGs_per_cluster <- split(lasc_annotation$OG, lasc_annotation$cluster)
lasc_OGs_per_cluster <- lapply(lasc_OGs_per_cluster, function(x) x[!is.na(x)])

# ----- clustering res.0.1 ------
# select columns from cluster biomarkers and ensure correct sorting
markers_res.0.1 <- read.csv("C:/Users/Cleiti/Desktop/Seurat_spiders/clavata/clavata_results/markers_res.0.1.csv", header = TRUE, sep = ",")
markers_res.0.1 <- markers_res.0.1[, c("p_val_adj", "cluster", "gene")]
markers_res.0.1 <- markers_res.0.1[order(markers_res.0.1$cluster, markers_res.0.1$p_val_adj), ] # making sure that it's sorted the right way
# subset top 50 genes per cluster
markers_res.0.1 <- as.data.frame(markers_res.0.1 %>% 
                                   group_by(cluster) %>%
                                   slice_min(order_by = p_val_adj, n = 50) %>%  
                                   ungroup())
# add OG for each gene using clavata annotation file (translation table)
markers_res.0.1$OG <- NA
for (i in 1:length(markers_res.0.1[,1])) {
  gene <- markers_res.0.1[i,3]
  if ((!is.na(gene)) && gene %in% OGs[,1]) {
    row <- which(OGs[,1] == gene)
    OG <- OGs[row, 2]
    markers_res.0.1[i,4] <- OG
  }
}
# create a nested list of OGs in each cluster
clavata_res.0.1_OGs_per_cluster <- split(markers_res.0.1$OG, markers_res.0.1$cluster)
clavata_res.0.1_OGs_per_cluster <- lapply(clavata_res.0.1_OGs_per_cluster, function(x) x[!is.na(x)])
# calculate overlap between each cluster in clavata res.0.1 and lasc
calculate_overlap <- function(genes1, genes2) {
  length(intersect(genes1, genes2))
}
overlap_matrix_res.0.1 <- outer(
  names(clavata_res.0.1_OGs_per_cluster), names(lasc_OGs_per_cluster), 
  Vectorize(function(cl1, cl2) calculate_overlap(clavata_res.0.1_OGs_per_cluster[[cl1]], lasc_OGs_per_cluster[[cl2]]))
)
# set rownames and colnames as cluster names
overlap_df_res.0.1 <- as.data.frame(overlap_matrix_res.0.1, row.names = names(clavata_res.0.1_OGs_per_cluster))
colnames(overlap_df_res.0.1) <- names(lasc_OGs_per_cluster)
# make a heatmap
pheatmap(overlap_df_res.0.1, 
         cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
         main = "Celltype overlaps: t.clavata res.0.1 vs l.sclopetarius",
         color = colorRampPalette(c("lightblue", "red"))(50))
