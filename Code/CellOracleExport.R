setwd("~/Documents/DPhil/Data/CHAPTER1/SINGLE_CELL/Chu_Cancer/CellOracle")

# Load libraries
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(slingshot)

# Load your data
seu <- readRDS("../Data/CD8_processed.rds")
sce <- readRDS("../Data/scesling.rds")

# Look at the pseudotime matrix
pseudotime_data <- slingPseudotime(sce)
seu$pseudotime_exhaustion <- pseudotime_data[Cells(seu), "Lineage2"]

FindVariableFeatures(seu)

# Convert to AnnData format
SaveH5Seurat(seu, filename = "tcell_data.h5Seurat")
Convert("tcell_data.h5Seurat", dest = "h5ad")
