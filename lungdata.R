# install.packages()
install.packages('tidyverse')
library(tidyverse)

install.packages('Seurat')
library(Seurat)

install.packages(c("Seurat", "tidyverse"))

install.packages('patchwork')
#Önceden işlenmiş scRNAseq ifade matrisini yükleyin

data <- readRDS("C:/Users/OZAN/Downloads/lung_ts.rds")

str(data)
head(data@meta.data)
View(data)

data

#Quality Check
VlnPlot(data, features = c("n_genes", "percent_mito", "n_counts"), ncol = 3, pt.size = 0.01)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(data, feature1 = "n_counts", feature2 = "percent_mito")
plot2 <- FeatureScatter(data, feature1 = "n_counts", feature2 = "n_genes")
plot1 + plot2

head(data@meta.data, 5)


plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#Normalizing the data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

#Değişken genleri belirleyin 
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
top10

#Verileri ölçekleme 
data <- ScaleData(data)

#Principal Component Analysis 
data <- RunPCA(data, features = VariableFeatures(object = data))
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)



data<- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:15)
ElbowPlot(data)

#Cluster cells by transcriptome
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
head(Idents(data), 5)

#Doğrusal olmayan boyut küçültmeyi çalıştırın (UMAP/tSNE)
data <- RunUMAP(data, dims = 1:10)
DimPlot(data, reduction = "umap")

saveRDS(data, file = "C:/Users/OZAN/Downloads/lung_ts.rds")

#deneme amaçlı clusterların markerları bulun
cluster0.markers <- FindMarkers(data, ident.1 = 2, min.pct = 0.25)
head(cluster0.markers, n = 5)

cluster24.markers <- FindMarkers(data, ident.1 = 2, min.pct = 0.25, logfc.threshold = 0.25)
head(cluster24.markers)
cluster24.markers

#Diferansiyel olarak ifade edilen özellikleri bulma (küme biyobelirteçleri)
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Her kümede belirli biyobelirteç genlerinin ifade dağılımı
VlnPlot(data, features = c("GNLY", "NKG7", "RGS1", "IL7R", "XCL2", "CD3D"))

FeaturePlot(data, features = c("GNLY", "NKG7", "RGS1", "IL7R", "XCL2", "CD3D"))


#Biyobelirteç gen ifadesinin ısı haritası görselleştirmesi
library(dplyr)
data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(data, features = top10$gene) + NoLegend()

#Kümelere hücre türü kimliği atama
new.cluster.ids <- c("Alveolar_Type1", "Alveolar_Type2", "Basal", "Blood_vessel", "B_cell_mature", "B_cell_naive", "Ciliated", "DC_1", "DC_2", "DC_activated", "DC_Monocyte_Dividing", "DC_plasmacytoid", "Fibroblast", "Lymph_vessel", "Macrophage_Dividing", "Macrophage_MARCOneg", "Macrophage_MARCOpos", "Mast_cells", "Monocyte", "Muscle_cells", "NK", "NK_Dividing", "Plasma_cells", "Secretory_club", "T_CD4", "T_CD8_CytT", "T_cells_Dividing", "T_regulatory")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

install.packages('rmarkdown')
rmarkdown::render("lungdata.R")

install.packages('tutorial')