# lungdataanalysis
scRNA-seq analysis 
### **Background**

The Human Cell Atlas is a large international collaborative effort to map all cell types of the human body. Single-cell RNA sequencing can generate high-quality data for the delivery of such an atlas. However, delays between fresh sample collection and processing may lead to poor data and difficulties in experimental design.

**Results**

This study assesses the effect of cold storage on fresh healthy spleen, esophagus, and lung from ≥ 5 donors over 72 h. We collect 240,000 high-quality single-cell transcriptomes with detailed cell type annotations and whole genome sequences of donors, enabling future eQTL studies. Our data provide a valuable resource for the study of these 3 organs and will allow cross-organ comparison of cell types.

We see little effect of cold ischemic time on cell yield, total number of reads per cell, and other quality control metrics in any of the tissues within the first 24 h. However, we observe a decrease in the proportions of lung T cells at 72 h, higher percentage of mitochondrial reads, and increased contamination by background ambient RNA reads in the 72-h samples in the spleen, which is cell type specific.

**Conclusions**

In conclusion, we present robust protocols for tissue preservation for up to 24 h prior to scRNA-seq analysis. This greatly facilitates the logistics of sample collection for Human Cell Atlas or clinical studies since it increases the time frames for sample processing.

**Good scRNA-seq data quality after cold storage**

The transplant surgeon assessed each organ as overall healthy in appearance. Whole genome sequencing (WGS) was carried out for each individual, confirming that none of the study participants displayed gross genomic abnormalities. 

.

This confirmed all tissue sections as healthy, except one donor with possible lung hypertension.

 **Annotation of cell types**

In the lung, **57,020 cells** passed quality control and represented 25 cell types. We detected `ciliated`, `alveolar types 1` and `2` cells, as well as `fibroblast`, `muscle`, and `endothelial` cells both from blood and lymph vessels. The cell types identified from the immune compartment included `NK`, `T`, and `B cells`, as well as two types of `macrophages`, `monocytes`, and `dendritic cells (DC)`. Multiple DC populations such as conventional `DC1`, `plasmacytoid DC (pcDC)`, and `activated DC` were detected and constituted 0.3% (163 cells), 0.08% (46 cells), and 0.2% (122 cells) of all cells, respectively. Lung club cell marker genes are detected in a small number of cells, but our clustering algorithm did not recognize these cells as a separate cluster . All donors contributed to every cluster. Dividing cells formed separate clusters for T cells, DC, monocytes, NK, and macrophages.

The gene expression count matrices from Cell Ranger output were used to perform sequential clustering of cells from either whole tissues or particular subclusters. The cell type identities of the clusters were determined and annotated by observation of expression of known cell type markers.

  **Data source**

https://www.tissuestabilitycellatlas.org/ (lung)

  **Paper**
  
Madissoon, E., Wilbrey-Clark, A., Miragaia, R.J. et al. scRNA-seq assessment of the human lung, spleen, and esophagus tissue stability after cold preservation. Genome Biol 21, 1 (2020). https://doi.org/10.1186/s13059-019-1906-x

