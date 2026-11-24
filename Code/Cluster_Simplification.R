library(tidyverse)
library(Seurat)


data <- readRDS('Data/CD8.rds')

Idents(data) <- 'cell.type'

DimPlot(data)

unique(data$cell.type)[!unique(data$cell.type) %in% names(cluster_mapping)]

cluster_mapping <- c(
    "CD8_c0_t-Teff"        = "Early Activation",
    "CD8_c1_Tex"          = "Exhausted",
    "CD8_c2_Teff"         = "Effector",
    "CD8_c3_Tn"           = "Naive",
    "CD8_c4_Tstr"         = "Stress Response",
    "CD8_c5_Tisg"         = "IFN-Stimulated",
    "CD8_c6_Tcm"          = "Central Memory",
    "CD8_c7_p-Tex"        = "Precursor Exhausted",
    "CD8_c8_Teff_KLRG1"   = "Effector",
    "CD8_c9_Tsen"         = "Senescent",
    "CD8_c10_Teff_CD244"  = "Effector",
    "CD8_c11_Teff_SEMA4A" = "Effector",
    "CD8_c12_Trm"         = "Resident Memory",
    "CD8_c13_Tn_TCF7"     = "Naive"
)

data$simplecluster <- unname(cluster_mapping[as.character(data$cell.type)])
table(data$simplecluster)
Idents(data) <- 'simplecluster'
DimPlot(data)

markers <- FindAllMarkers(data, logfc.threshold = 0.5, min.diff.pct = 0.2, only.pos = T, return.thresh = 0.05)

exhaustion <- c('HAVCR2','PDCD1','LAG3','CTLA4','TIGIT','TOX','ENTPD1','LAYN')
effector <- c('GZMB','GZMH','GNLY','PRF1','IFNG','CX3CR1','NKG7')
naive <- c('CCR7','SELL','LEF1','TCF7','IL7R','CD28','CD27')
preexhausted <- c('TCF7','XCL1','XCL2','SLAMF6','CXCR5')
memory <- c('GZMK','IL7R','EOMES','CD44','CCR5','CXCR3')
rm <- c('ITGAE','CD69','CXCR6','ITGA1','ZNF683')
proliferation <- c('MKI67','TOP2A','PCNA','CDK1')

allmarkers <- c(exhaustion,effector,naive,preexhausted,memory,rm,proliferation)

markerDEGs <- markers %>%
    group_by(cluster) %>%
    filter (gene %in% allmarkers, avg_log2FC > 0.5, p_val < 0.1, pct.1 > 0.5)
