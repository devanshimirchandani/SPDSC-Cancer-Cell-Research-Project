library(Seurat)
library(dplyr)
library(ggplot2)


gene_data <- CrC_counts
metadata <- CrC_metadata

# Check orig.ident
unique(CrC_metadata$orig.ident)

crc_seurat[["percent.mt"]] <- PercentageFeatureSet(crc_seurat, pattern = "^MT-")

CrC_metadata_subset <- subset(CrC_metadata, orig.ident %in% c("B001-A-301", "F007", "CRC1_8810"))

# Subset the counts data based on the cells in the filtered metadata
CrC_counts_subset <- CrC_counts[, rownames(CrC_metadata_subset)]

# Create the Seurat object with the subset data
crc_seurat <- CreateSeuratObject(counts = CrC_counts_subset, meta.data = CrC_metadata_subset, project = "CRC_analysis", min.cells = 3, min.features = 200)
crc_seurat


VlnPlot(crc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

crc_seurat <- subset(crc_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)

# Normalisation
crc_seurat <- NormalizeData(crc_seurat)

# Feature Scatter plots
plot1 <- FeatureScatter(crc_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(crc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

crc_seurat <- FindVariableFeatures(crc_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(crc_seurat), 10)

# Plot variable features with and without labels
plot_a <- VariableFeaturePlot(crc_seurat)
plot_b <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot_a + plot_b

# Scaling data
all.genes <- rownames(crc_seurat)
crc_seurat <- ScaleData(crc_seurat, features = all.genes)

## PCA
crc_seurat <- RunPCA(crc_seurat, features = VariableFeatures(crc_seurat))
VizDimLoadings(crc_seurat, dims = 1:2, reduction = "pca") 
DimPlot(crc_seurat, reduction = "pca")

## UMAP
crc_seurat <- RunUMAP(crc_seurat, dims = 1:10) 
DimPlot(crc_seurat, reduction = "umap")


