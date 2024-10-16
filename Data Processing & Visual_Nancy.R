library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load data
CrC_counts <- readRDS("CrC_counts.Rds")
CrC_metadata <- readRDS("CrC_metadata.Rds")

# Create Seurat object
crc_seurat <- CreateSeuratObject(counts = CrC_counts)
crc_seurat <- AddMetaData(object = crc_seurat, metadata = as.data.frame(CrC_metadata))

# Check donors and disease states
print(table(crc_seurat$Donor, crc_seurat$DiseaseState))

# Function to process Seurat object
process_seurat <- function(seurat_obj) {
  seurat_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)
}

# Process the entire dataset
crc_seurat <- process_seurat(crc_seurat)

# Split the dataset based on disease state, ensuring all donors are included
Idents(crc_seurat) <- crc_seurat$DiseaseState
seurat_list <- SplitObject(crc_seurat, split.by = "DiseaseState")

# Check donors in each split
for(name in names(seurat_list)) {
  print(paste("Donors in", name, "group:"))
  print(table(seurat_list[[name]]$Donor))
}

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                  anchor.features = 2000, 
                                  reduction = "rpca")

# Integrate data
integrated_seurat <- IntegrateData(anchorset = anchors, new.assay.name = "integrated")

# Check donors after integration
print(table(integrated_seurat$Donor))

# Set the default assay to "integrated" and scale the data
DefaultAssay(integrated_seurat) <- "integrated"
integrated_seurat <- ScaleData(integrated_seurat, verbose = FALSE)

# Perform PCA, find neighbors, clustering, and run UMAP on the integrated data
integrated_seurat <- integrated_seurat %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:30)

# PCA plot colored by donor labels
DimPlot(integrated_seurat, reduction = "pca", group.by = "Donor")

# UMAP plot colored by donor labels
DimPlot(integrated_seurat, reduction = "umap", group.by = "Donor")

# UMAP plot colored by cluster
DimPlot(integrated_seurat, reduction = "umap", group.by = "seurat_clusters")

# UMAP plot split by DifferentialGroup
DimPlot(integrated_seurat, reduction = "umap", split.by = "DifferentialGroup")

# UMAP plot split by DiseaseState
DimPlot(integrated_seurat, reduction = "umap", split.by = "DiseaseState")



# Visualize UMAP, colored by cell type (e.g., "seurat_annotations") and condition ("DiseaseState")
DimPlot(integrated_seurat, reduction = "umap", group.by = c("DiseaseState", "seurat_clusters"))

# Side-by-Side Visualization of Disease States
DimPlot(integrated_seurat, reduction = "umap", split.by = "DiseaseState")
