suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(ComplexHeatmap)
    library(ggplot2)
    library(gridExtra)
    library(Hmisc)
    library(STRINGdb)
    library(igraph)
    library(SingleR)
    library(ggpubr)
    library(ggrepel)
})
# Feature lists ----
AllIRs <- as.character(read.delim('../../IRListGeneNames.txt')[,1])
IRs <- c('TIGIT','HAVCR2','PDCD1','CTLA4','BTLA','LAG3')


# load data ----

data<-readRDS('Data/CD8.rds')
data
colnames(data@meta.data)
table(data@meta.data$cell.type)
Idents(data)<-'cell.type'
dir.create('Plots/')

Features(data)
Idents(data)
VlnPlot(data, features = 'HAVCR2', group.by = 'cell.type')
DimPlot(data, group.by = 'cell.type')

DimPlot(data, group.by = 'TissueType')

table(Idents(data))
data@meta.data
data@meta.data$cell.type <- gsub("CD8_", "", data@meta.data$cell.type)
Idents(data) <- 'cell.type'
Idents(data)
data@meta.data$cell.type


# Pseudobulk + PCA ----

library(Seurat)
library(dplyr)
library(ggplot2)



# 1. Create pseudobulk expression matrix
# Aggregate expression by cell type
pseudobulk_matrix <- AverageExpression(
    data, 
    group.by = "cell.type",
    return.seurat = FALSE
)$RNA

# 2. Filter out genes with zero or constant expression
# Remove genes with zero expression across all cell types
genes_to_keep <- rowSums(pseudobulk_matrix) > 0

# Remove genes with constant expression (zero variance)
gene_variance <- apply(pseudobulk_matrix, 1, var)
genes_with_variance <- gene_variance > 0

# Combine filters
final_genes_to_keep <- genes_to_keep & genes_with_variance

# Filter the matrix
pseudobulk_filtered <- pseudobulk_matrix[final_genes_to_keep, ]

print(paste("Original genes:", nrow(pseudobulk_matrix)))
print(paste("Genes after filtering:", nrow(pseudobulk_filtered)))

# 3. Log-transform the filtered pseudobulk data
pseudobulk_log <- log2(pseudobulk_filtered + 1)

# Transpose so cell types are rows and genes are columns
pseudobulk_t <- t(pseudobulk_log)

# 3. Perform PCA
pca_result <- prcomp(pseudobulk_t, scale = TRUE, center = TRUE)

# 4. Create PCA plot
# Extract PC coordinates
pca_coords <- data.frame(
    cell_type = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
)

# Calculate variance explained
variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Create the plot
ggplot(pca_coords, aes(x = PC1, y = PC2, color = cell_type, label = cell_type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(vjust = -0.5, hjust = 0.5, size = 3, colour = 'black') +
    labs(
        x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
        title = "PCA of Pseudobulk Expression by Cell Type"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)
    )

ggsave('Plots/Pseudobulk/PCAplot.pdf')



# Feature overlay plots ----
dir.create('Plots/Dimplots/')

DimPlot(data, group.by = 'cell.type')
ggsave('Plots/dimplotbycluster.pdf')
DimPlot(data, group.by = 'CancerType')
ggsave('Plots/dimplotbycancer.pdf')
DimPlot(data, group.by = 'TissueType')
ggsave('Plots/dimplotbytissue.pdf')

gene <- 'HAVCR2'
FeaturePlot(data, features = gene, blend.threshold = 0.5, label = T, max.cutoff = 5, label.size = 3, label.color = 'black', repel = T) 
ggsave(paste0('Plots/Dimplots/',gene,'.pdf'))

FeaturePlot(data, features = IRs ,blend.threshold = 0.5, min.cutoff = 0, max.cutoff = 5, label = T, label.size = 3, label.color = 'black', repel = T) 
ggsave('Plots/Dimplots/allIRs.pdf', width =10, height = 13)

for (i in (1:length(IRs))){
    gene <- IRs[i]
    p0<-DimPlot(data, label = T, repel = T) + ggtitle('T Cell Clusters') + theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    p1<-FeaturePlot(data, features = gene, min.cutoff = 0, max.cutoff = 5, label = F, repel = T) + ggtitle(paste0(gene, ' expression')) + theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    p<-p0|p1
    ggsave(paste0('Plots/Dimplots/',gene,'_dimplot.pdf'), width = 15)
}

# Violin plots ----
dir.create('Plots/Vlnplots')


excluster<-unique(grep('c1_Tex',data@meta.data$cell.type, value=T))
pexcluster <- unique(grep('p-Tex',data@meta.data$cell.type, value=T))
naivecluster <- unique(grep('Tn',data@meta.data$cell.type, value=T))
gene <- 'CD8A'
gene_data <- FetchData(data, vars = c(gene, "TissueType"))
colnames(gene_data) <- c("expression", "cluster")
# gene_data <- gene_data[gene_data$cluster != 'Healthy donor',]
# gene_data$variable<- 'All Other Cells'
# gene_data$variable[gene_data$cluster %in% excluster] <- 'Terminally Exhausted Cells'
# gene_data$variable[gene_data$cluster %in% pexcluster] <- 'Precursor Exhausted Cells'
# gene_data$variable[gene_data$cluster %in% naivecluster] <- 'Naive Cells'



# Create the plot
ggplot(gene_data, aes(x = factor(cluster), y = expression)) +
    geom_violin(aes(fill = factor(cluster)), alpha = 0.7, bounds = c(0,5)) +
    geom_boxplot(width = 0.2, alpha = 0.5, outliers = F) +
    labs(title = paste0(gene, " Expression by Cell Type"),
         x = "Cluster", y = "Expression Level") +
    theme_bw() +
    theme(legend.position = "none") 
ggsave(paste0('Plots/Vlnplots/',gene,'_vlnplot_by_celltype.pdf'), width =20)





for (i in (1:length(IRs))){
    gene <- IRs[i]
    pdf(paste0('Plots/Vlnplots/', gene, 'vlnplot.pdf'), width = 10)
    p<-VlnPlot(data, features = gene, pt.size = 0.2) + theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5))
    print(p)
    dev.off()    
}

# Heatmap of feature lists ----
dir.create('Plots/Heatmaps')


plottingdata <- subset(data, downsample = 100)
gene_matrix <- GetAssayData(plottingdata, assay = "RNA", slot = "data")
gene_matrix <- as.matrix(gene_matrix[IRs, ])
cell_annotations <- factor(plottingdata$cell.type)

# Plot heatmap with hierarchical clustering 

pdf('Plots/Heatmaps/IRs.pdf')
p<- ComplexHeatmap::Heatmap(gene_matrix, 
            name = 'IR Expression',
            cluster_columns = FALSE,
            show_column_dend = FALSE,
            cluster_rows = FALSE,
            show_column_names = FALSE,  
            column_split = cell_annotations,
            column_title_gp = gpar(fontsize = 8),
            row_names_gp = gpar(fontsize = 6),
)
draw(p)
dev.off()

# Dotplots ----
dir.create('Plots/Dotplots')

DotPlot(data, features = IRs,col.min = 0, scale.by = 'radius', group.by = 'cell.type', scale = T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Plots/Dotplots/IRdotplot.pdf')











# HAVCR2 Coexpression ----


library(Seurat)
library(Matrix)
library(dplyr)

# Step 1: Filter cells by HAVCR2 expression first (keep top 10%)
print("Step 1: Filtering cells by HAVCR2 expression...")

# Check if HAVCR2 is in your dataset
if(!"HAVCR2" %in% rownames(data)) {
    stop("HAVCR2 not found in the dataset. Check gene name spelling.")
}

# Extract HAVCR2 expression from data layer
havcr2_expr <- data@assays$RNA$data["HAVCR2", ]

# Calculate 90th percentile threshold
havcr2_threshold <- quantile(havcr2_expr, 0.9)
print(paste("HAVCR2 90th percentile threshold:", round(havcr2_threshold, 3)))

# Identify cells in top 10% of HAVCR2 expression
high_havcr2_cells <- names(havcr2_expr)[havcr2_expr > havcr2_threshold]

print(paste("Original cells:", length(havcr2_expr)))
print(paste("Cells after HAVCR2 filtering:", length(high_havcr2_cells)))

# Create cell-filtered Seurat object first
print("Creating cell-filtered dataset...")
data_cell_filtered <- subset(data, cells = high_havcr2_cells)

print(paste("Dataset after cell filtering:", ncol(data_cell_filtered), "cells x", nrow(data_cell_filtered), "genes"))

# Step 2: Now remove low variability genes (much smaller dataset now)
print("Step 2: Filtering low variability genes...")

# Extract the data layer from the cell-filtered object
expr_matrix <- data_cell_filtered@assays$RNA$data

# Calculate variance for each gene (now much more manageable)
print("Calculating gene variances...")
gene_vars <- apply(expr_matrix, 1, var)
names(gene_vars) <- rownames(expr_matrix)

# Remove genes with zero variance
genes_with_variance <- gene_vars > 0

# Keep top variable genes (top 50% most variable)
var_threshold <- quantile(gene_vars[genes_with_variance], 0.5)
high_var_genes <- names(gene_vars)[gene_vars > var_threshold]

print(paste("Original genes:", length(gene_vars)))
print(paste("Genes after variance filtering:", length(high_var_genes)))

# Step 3: Create final filtered Seurat object
print("Creating final filtered dataset...")

# Subset by genes (cells already filtered)
data_filtered <- subset(data_cell_filtered, features = high_var_genes)

print(paste("Final dataset size:", ncol(data_filtered), "cells x", nrow(data_filtered), "genes"))

# Optional: Check cell type composition after filtering
if("cell.type" %in% names(data_filtered@meta.data)) {
    print("Cell type composition after filtering:")
    print(table(data_filtered$cell.type))
}

# Step 4: Run coexpression analysis
print("Running coexpression analysis...")

# Extract expression matrix from data layer for correlation analysis
expr_for_cor <- as.matrix(data_filtered@assays$RNA$data)

# Transpose so genes are columns (required for correlation)
expr_for_cor_t <- t(expr_for_cor)

# Calculate correlations with HAVCR2
print("Calculating correlations with HAVCR2...")
havcr2_correlations <- cor(expr_for_cor_t[, "HAVCR2"], expr_for_cor_t, method = "pearson")

# Convert to data frame and sort
cor_results <- data.frame(
    gene = colnames(havcr2_correlations),
    correlation = as.vector(havcr2_correlations),
    stringsAsFactors = FALSE
) %>%
    arrange(desc(abs(correlation))) %>%
    filter(gene != "HAVCR2")  # Remove HAVCR2 itself

# Display top results
print("Top 20 positively correlated genes with HAVCR2:")
print(head(cor_results %>% filter(correlation > 0), 20))

print("Top 20 negatively correlated genes with HAVCR2:")
print(head(cor_results %>% filter(correlation < 0), 20))

# Plot correlation distribution
ggplot(cor_results, aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.9) +
    labs(title = "Distribution of HAVCR2 Gene Correlations",
         x = "Pearson Correlation with HAVCR2",
         y = "Number of Genes") +
    theme_bw()

ggsave('Plots/Coexpression/HAVCR2_Distribution.pdf')

# Plot top correlated genes
top_genes <- rbind(
    head(cor_results %>% filter(correlation > 0), 20)
)

ggplot(top_genes, aes(x = reorder(gene, correlation), y = correlation)) +
    geom_col(alpha = 0.9, fill = 'skyblue', color = 'black') +
    coord_flip() +
    labs(title = "Top 20 Genes Correlated with HAVCR2",
         x = "Gene",
         y = "Pearson Correlation") +
    theme_bw()

ggsave('Plots/Coexpression/HAVCR2_PositiveCorrelations.pdf')

# Save results
write.csv(cor_results, "HAVCR2_coexpression_results.csv", row.names = FALSE)
print("Results saved to HAVCR2_coexpression_results.csv")

