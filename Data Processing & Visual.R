library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

gene_data <- CrC_counts
metadata <- CrC_metadata

unique(CrC_metadata$Donor)

crc_seurat <- CreateSeuratObject(counts = CrC_counts)
crc_seurat <- AddMetaData(object = crc_seurat, metadata = as.data.frame(CrC_metadata))

# View the basic structure of the Seurat object
crc_seurat

length(Cells(crc_seurat))
length(crc_seurat@meta.data$Donor)

table(crc_seurat@meta.data$Donor)

cells_by_donor <- split(Cells(crc_seurat), crc_seurat@meta.data$Donor)

set.seed(123)  # Set seed for reproducibility 
sampled_cells <- unlist(lapply(cells_by_donor, function(cells) { sample(cells, size = min(5000, length(cells)), replace = FALSE)}))

sampled_crc_seurat <- subset(crc_seurat, cells = sampled_cells)
sampled_crc_seurat

sampled_crc_seurat[["percent.mt"]] <- PercentageFeatureSet(sampled_crc_seurat, pattern = "^MT-")

VlnPlot(sampled_crc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

sampled_crc_seurat <- subset(sampled_crc_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)

sampled_crc_seurat <- NormalizeData(sampled_crc_seurat)

plot1 <- FeatureScatter(sampled_crc_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Donor")
plot2 <- FeatureScatter(sampled_crc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Donor")

plot1
plot2


sampled_crc_seurat <- FindVariableFeatures(sampled_crc_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sampled_crc_seurat), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(sampled_crc_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1
plot2

all.genes <- rownames(sampled_crc_seurat)
sampled_crc_seurat <- ScaleData(sampled_crc_seurat, features = all.genes)

sampled_crc_seurat

sampled_crc_seurat <- RunPCA(sampled_crc_seurat, features = VariableFeatures(sampled_crc_seurat))

VizDimLoadings(sampled_crc_seurat, dims = 1:2, reduction = "pca") 

DimPlot(sampled_crc_seurat, reduction = "pca", group.by = "Donor")

sampled_crc_seurat <- RunUMAP(sampled_crc_seurat, dims = 1:10) 

DimPlot(sampled_crc_seurat, reduction = "umap", group.by = "Donor")



