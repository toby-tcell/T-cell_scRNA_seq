library(Seurat)
library(tidyverse)
library(ggrepel)
library(slingshot)
library(ggpubr)
library(colorspace)
library(ggnewscale)
library(RColorBrewer)


setwd("~/Documents/DPhil/Data/CHAPTER3/SINGLE_CELL/Chu_Cancer")

plot_path <- 'General_Plots/'
dir.create(plot_path, showWarnings = F)

data <- readRDS('Data/CD8.rds')
data$cell.type <- gsub('CD8_','', data$cell.type)
colnames(data@meta.data)

table(data$TissueType)

tissue_mapping <- c(
    "Healthy donor"        = "Healthy Donor",
    "Metastatic tumor tissue"          = "Metastatic Tumour",
    "Primary tumor tissue"         = "Primary Tumour",
    "Uninvolved normal tissue"           = "Adjacent Tissue"
)

data$TissueType <- unname(tissue_mapping[as.character(data$TissueType)])
table(data$TissueType)
data$TissueType <- factor(data$TissueType, levels = c('Healthy Donor','Adjacent Tissue','Primary Tumour','Metastatic Tumour'))


# Cell filtering 

data[['percent.mt']] <- PercentageFeatureSet(data, assay = 'RNA', pattern = "^MT-")
data[['percent.rb']] <- PercentageFeatureSet(data, assay = 'RNA', pattern = "^RP[SL]")
data[['percent.hb']] <- PercentageFeatureSet(data, assay = 'RNA', pattern = "^HB[AB]")
data$nCount_RNA <- Matrix::colSums(data@assays$RNA$counts)
data$nFeature_RNA <- Matrix::colSums(data@assays$RNA$counts > 0)

head(data@meta.data)

qc.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb")

for (i in 1:length(qc.features)){
    feature = qc.features[i]
    VlnPlot(data, features = feature, group.by = 'cell.type') +
        labs(y = NULL) +
        theme(legend.position = 'none')
    ggsave(paste0(plot_path, feature, '.pdf'), width = 5, height = 5)
    
    
}


data_qc <- subset(data, 
                  subset = nFeature_RNA >200 &
                      nFeature_RNA < 6000 &
                      percent.mt < 5 & 
                      percent.hb < 2
                  )

data <- data_qc

# UMAP plot of existing UMAP
UMAPPlot(data, group.by = 'cell.type', label= T, label.size = 3.5, repel = T, stroke.size =0.25, raster = F, label.box = T) + 
    ggtitle(label = NULL) +
    theme(legend.position = 'none')
ggsave(paste0(plot_path, 'UMAP_labels_original.pdf'), width = 6, height = 6)

UMAPPlot(data, group.by = 'cell.type', label= F,  stroke.size =0.5, raster = T) + 
    ggtitle(label = NULL) +
    theme(legend.position = 'none')
ggsave(paste0(plot_path, 'UMAP_fig3.0.pdf'), width = 6, height = 6)


# Cluster simplification

cluster_mapping <- c(
    "c0_t-Teff"        = "Early Activation",
    "c1_Tex"          = "Exhausted",
    "c2_Teff"         = "Effector",
    "c3_Tn"           = "Naive",
    "c4_Tstr"         = "Stress Response",
    "c5_Tisg"         = "IFN-Stimulated",
    "c6_Tcm"          = "Central Memory",
    "c7_p-Tex"        = "Precursor Exhausted",
    "c8_Teff_KLRG1"   = "Effector",
    "c9_Tsen"         = "Senescent", 
    "c10_Teff_CD244"  = "Effector",
    "c11_Teff_SEMA4A" = "Effector",
    "c12_Trm"         = "Resident Memory",
    "c13_Tn_TCF7"     = "Naive"
)

data$simplecluster <- unname(cluster_mapping[as.character(data$cell.type)])
table(data$simplecluster)
Idents(data) <- 'simplecluster'



# UMAP new cell types

UMAPPlot(data, group.by = 'simplecluster', label= T, label.size = 4, repel = T, stroke.size =0.25, raster = F, label.box = T) + 
    ggtitle(label = NULL) +
    theme(legend.position = 'none')
ggsave(paste0(plot_path, 'UMAP_labels_new.pdf'), width = 6, height = 6)





# CELL TYPE DISTRIBUTIONS

Celltypedistributions <- data.frame(data@meta.data) %>%
    group_by(across(all_of(c('simplecluster', 'TissueType')))) %>%
    summarise(n = n(), .groups = 'drop') %>%
    group_by(across(all_of('TissueType'))) %>%
    mutate(pct.cell = prop.table(n)*100)


ggplot(Celltypedistributions, aes(x = TissueType, y = pct.cell, fill = simplecluster)) +
    geom_bar(stat = 'identity', show.legend = F, color = 'black', linewidth = 0.2) +
    ylab('%Subset') +
    geom_text(aes(label = ifelse(pct.cell < 3, '' , paste0(simplecluster))), 
              position = position_stack(vjust = 0.5), size = 2) +
    xlab(NULL) +
    scale_fill_discrete_diverging(palette = 'Blue-Red') +
    scale_x_discrete(guide = guide_axis(angle = 0)) + 
    theme_bw()
ggsave(paste0(plot_path, 'Celltypebytissuetype.pdf'), height = 5, width = 5)

# Markers DotPlot

markergenes <- c('TCF7','SELL','IL7R','CD69','NKG7','GNLY','FGFBP2','CCR7','CD44','HSPA1A','HSPA1B','IFIT1','MX1','SLAMF6','PDCD1','CTLA4','LAG3','KLRC4')
idents <- c('Naive','Early Activation','Effector','Resident Memory','Central Memory','Stress Response','IFN-Stimulated','Precursor Exhausted','Exhausted','Senescent')
Idents(data) <- factor(Idents(data), levels = idents)
table(Idents(data))
DotPlot(data, features = markergenes) + scale_x_discrete(guide = guide_axis(angle = 45)) + ylab(NULL) + xlab(NULL)
ggsave(file.path(plot_path, 'Marker_DotPlot.pdf'), height = 7, width = 12)


# UMAP Plots of key markers

FeaturePlot(data, features = 'SELL', raster = F)


# Functional scores per cell

all_signatures <- list(
    Cytotoxicity = c("PRF1", "GZMB", "GZMA", "GZMK", "GNLY", "NKG7", "FASLG", "GZMH",'GZMM','FAS','TNFSF10','LAMP1'),
    Activation = c("CD69", "IL2RA", "CD44", "ICOS", "TNFRSF9", "TNFRSF4", "KLRG1",'CD38','HLA-DR','FOS','JUN','NFATC1'),
    Exhaustion = c("PDCD1", "LAG3", "TIGIT","LAYN",'CTLA4','BTLA','CD244','CD160'),
    Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "CCNB1", "CDK1", "CCNA2"),
    Survival = c("BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "BBC3", "BCL2L11"),
    Memory_Formation = c("TCF7", "LEF1", "CCR7", "IL7R", "CD62L",'BACH2','ID3','FOXO1'),
    Tissue_Residency = c("CD69", "ITGAE", "CXCR6", "CD101",'ZNF683','ITGA1','RGS1','RGS2','RUNX3'),
    Proteotoxic_Stress = c("HSPA1A", "HSP90AA1", "DDIT3", "ATF4", "XBP1", "HSPA5", "NFE2L2", "SQSTM1"),
    Glycolysis = c("SLC2A1", "HK2", "PFKP", "ALDOA", "TPI1", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA"),
    OXPHOS = c("NDUFA1", "NDUFS1", "SDHA", "SDHB", "UQCRC1", "COX4I1", "ATP5F1A", "ATP5F1B", "MT-CO1", "MT-ND1"),
    Fatty_Acid_Oxidation = c("CPT1A", "ACADVL", "HADHA", "ACADM", "ACOX1", "EHHADH", "CPT2"),
    Fatty_Acid_Synthesis = c("ACLY", "ACACA", "FASN", "SCD1", "ELOVL6", "ACSS2"),
    Glutamine_Metabolism = c("SLC1A5", "GLS", "GLUD1", "GOT1", "GPT2", "ASNS", "GLS2"),
    Serine_Glycine = c("PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD2"),
    Pentose_Phosphate = c("G6PD", "PGD", "RPIA", "RPE", "TALDO1", "TKT", "RBKS"),
    TCA_Cycle = c("CS", "ACO2", "IDH2", "OGDH", "SUCLA2", "SUCLG1", "SDHA", "FH", "MDH2"),
    Amino_Acid_Transport = c("SLC1A5", "SLC7A5", "SLC3A2", "SLC38A2", "SLC6A14", "SLC7A1"),
    Mitochondrial_Biogenesis = c("TFAM", "NRF1", "PPARGC1A", "PPARGC1B", "ESRRA", "GABPA"),
    mTOR_Signaling = c("MTOR", "RPS6KB1", "EIF4EBP1", "AKT1", "TSC1", "TSC2", "RHEB"),
    AMPK_Signaling = c("PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "STK11")
)

for (signature_name in names(all_signatures)) {
    # Check which genes are present in the dataset
    genes_present <- all_signatures[[signature_name]][all_signatures[[signature_name]] %in% rownames(data)]
    
    if (length(genes_present) >= 3) {  # Need at least 3 genes for meaningful score
        data <- AddModuleScore(
            object = data,
            features = list(genes_present),
            name = paste0(signature_name, "_Score"),
            assay = "RNA"
        )
        print(paste("Calculated score for", signature_name, "using", length(genes_present), "genes"))
    } else {
        print(paste("Skipping", signature_name, "- insufficient genes present (", length(genes_present), ")"))
    }
}

analysis_df <- data.frame(
    cell_id = colnames(data),
    cell_type = data@meta.data$simplecluster,
    stringsAsFactors = FALSE
)

for (signature_name in names(all_signatures)) {
    score_col <- paste0(signature_name, "_Score1")  # Seurat adds "1" suffix
    if (score_col %in% colnames(data@meta.data)) {
        analysis_df[[paste0(signature_name, "_Score")]] <- data@meta.data[[score_col]]
    }
}

idents <- c('Naive','Early Activation','Effector','Resident Memory','Central Memory','Stress Response','IFN-Stimulated','Precursor Exhausted','Exhausted','Senescent')
analysis_df$cell_type <- factor(analysis_df$cell_type, levels = idents)

module_plot_path <- 'General_Plots/Modules/'
dir.create(module_plot_path, showWarnings = F)

# Create boxplots by cell type
for (signature_name in names(all_signatures)) {
    score_col <- paste0(signature_name, "_Score")
    if (score_col %in% colnames(analysis_df)) {
        p_box <- ggplot(analysis_df, aes_string(x = 'cell_type', y = score_col)) +
            geom_boxplot(aes(fill = cell_type), alpha = 0.8, outliers = F,) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none") +
            labs(x = NULL, y = paste(signature_name, "Score")) 
        ggsave(paste0(module_plot_path, tolower(signature_name), "_by_celltype.pdf"), height = 4, width = 5, p_box)
    }
}

# Check 


# Save data into RDS file ----

colnames(data@meta.data)
saveRDS(data, 'Data/CD8_processed.rds')









