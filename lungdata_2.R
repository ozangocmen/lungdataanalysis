# install.packages()
install.packages('tidyverse')
library(tidyverse)

install.packages('Seurat')
library(Seurat)

#or alternative installing 
install.packages(c("Seurat", "tidyverse"))
#the goal of that package, to combine separate ggplots into the same graphic.
install.packages('patchwork')
library(patchwork)

#load pre-processed dataset, my dataset name is data 
#(suggestion: use another variable name)
data <- readRDS("C:/Users/OZAN/Downloads/lung_ts.rds")

str(data)
head(data@meta.data)
View(data)

data
load
#our file is .rds format and as a Seurat object 

# 1.Quality Check and selecting cells 
percent_mito <- grep(pattern = "^MT-", x = rownames(x = data@data), value = TRUE)
data[["percent_mito"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("n_genes", "percent_mito", "n_counts"), ncol = 3, pt.size = 0.01)
head(data@meta.data, 5)
# feature comparison
plot1 <- FeatureScatter(data, feature1 = "n_counts", feature2 = "percent_mito")
plot2 <- FeatureScatter(data, feature1 = "n_counts", feature2 = "n_genes")
plot1 + plot2


# 2.Normalizing the data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

# 3.Feature selection 
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 4.Scaling the data 
data <- ScaleData(data)

# 5.Perform linear dimensional reduction
data <- RunPCA(data, features = VariableFeatures(object = data))
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

# 6.Determine the ‘dimensionality’ of the dataset
ElbowPlot(data)

# 7.Cluster the cells 
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
head(Idents(data), 5)

# 8.Run non-linear dimensional reduction (UMAP/tSNE)
data <- RunUMAP(data, dims = 1:10)
DimPlot(data, reduction = "umap")



cluster0.markers <- FindMarkers(data, ident.1 = 2, min.pct = 0.25)
head(cluster0.markers, n = 5)

cluster24.markers <- FindMarkers(data, ident.1 = 2, min.pct = 0.25, logfc.threshold = 0.25)
head(cluster24.markers)
cluster24.markers

# 9.Cluster biomarkers 
#find markers for every cluster compared to all remaining cells, report only the positive ones
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

data.markers

write_csv(data.markers, "C:\\Users\\OZAN\\Desktop\\bioinformatics\\Lung cancer CAR T gene therapy\\analysis\\1\\data.markers.csv")

library(data.table)

fwrite(data.markers, "C:\\Users\\OZAN\\Desktop\\bioinformatics\\Lung cancer CAR T gene therapy\\analysis\\1\\data.markers.csv")
#VlnPlot, shows expression probability distributions across clusters, and FeaturePlot
VlnPlot(data, features = c("GNLY", "NKG7", "RGS1", "IL7R", "XCL2", "CD3D"))

FeaturePlot(data, features = c("GNLY", "NKG7", "RGS1", "IL7R", "XCL2", "CD3D"))


#DoHeatmap, generates an expression heatmap for given cells and features.
library(dplyr)
data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(data, features = top10$gene) + NoLegend()

# 10.Assigning cell type identity to clusters
new.cluster.ids <- c("Alveolar_Type1", "Alveolar_Type2", "Basal", "Blood_vessel", "B_cell_mature", "B_cell_naive", "Ciliated", "DC_1", "DC_2", "DC_activated", "DC_Monocyte_Dividing", "DC_plasmacytoid", "Fibroblast", "Lymph_vessel", "Macrophage_Dividing", "Macrophage_MARCOneg", "Macrophage_MARCOpos", "Mast_cells", "Monocyte", "Muscle_cells", "NK", "NK_Dividing", "Plasma_cells", "Secretory_club", "T_CD4", "T_CD8_CytT", "T_cells_Dividing", "T_regulatory")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

sessionInfo()

save.image()
