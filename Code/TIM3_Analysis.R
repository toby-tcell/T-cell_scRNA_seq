
# Libraries ----
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(Matrix)
library(ggpubr)
library(ggrepel)
library(ggforce)
library(edgeR)
library(limma)
library(rstatix)



# Variables and data ----

gene <- 'HAVCR2'
output_path <- paste0('TIM3_plots/')
data_path <- 'Data/CD8_processed.rds'
pseudotime_path <- 'Data/scesling.rds'
dir.create(output_path, showWarnings = F)
cell_cluster_name <- 'simplecluster'

# Load data 

data <- readRDS(data_path)
Idents(data) <- cell_cluster_name

# Add exhaustion grouping

data$ExhaustionGroup <- ifelse(data$simplecluster == 'Exhausted', 'Exhausted', ifelse(data$simplecluster == 'Precursor Exhausted', 'Precursor Exhausted','Non-Exhausted'))

table(data$ExhaustionGroup)

# Making exploratory graph outputs ----

# Violin Plots

gene_data <- FetchData(data, vars = c(gene, cell_cluster_name, 'TissueType', 'ExhaustionGroup'))
colnames(gene_data) <- c("expression", "cluster",'tissue','exhaustion')

gene_data$exhaustion <- factor(gene_data$exhaustion, levels = c('Non-Exhausted','Precursor Exhausted','Exhausted'))
table(gene_data$exhaustion)

ggplot(gene_data, aes(x = exhaustion, y = expression)) +
    geom_violin(aes(fill = exhaustion), alpha = 1, trim = T, scale = 'width', width = 0.7) +
    geom_boxplot(width = 0.4, alpha = 0.5, outliers = F) +
    labs(y = "Normalized HAVCR2 Expression", x=NULL) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_discrete(guide = guide_axis(angle = 0)) 
ggsave(paste0(output_path, 'vlnplot_byexhaustion.pdf'), width =5, height = 5)



ggplot(gene_data, aes(x = tissue, y = expression)) +
    geom_violin(alpha = 0.5, scale = 'width', aes(fill = tissue), width = 0.6) +
    geom_boxplot(width = 0.4, alpha = 0.2, outliers = F) +
    labs(x = 'Tissue', y = "Normalized Expression") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    labs(x = NULL) 
ggsave(paste0(output_path, 'vlnplot_bytissue.pdf'), width =5, height = 5)


FeaturePlot(data, features = gene, label = F, raster = F) + ggtitle(NULL)
ggsave(paste0(output_path, 'featureoverlayplot.pdf'), height = 5, width = 5)


# Pseudobulk analysis ----



# Q1. Is HAVCR2 expressed higher in pre-exhausted and exhausted cells? ----

data$pb_id <- interaction(
    data$batch,
    data$CancerType,
    data$ExhaustionGroup,
    drop = T,
    sep = '--'
)
table(data$pb_id)
length(unique(data$pb_id))

pb_counts <- AggregateExpression(
    data,
    assays = 'RNA',
    slot = 'counts',
    group.by = 'pb_id',
    return.seurat = F
)$RNA %>%
    round()


pb_meta <- data@meta.data %>%
    as.data.frame() %>%
    distinct(pb_id, batch, TissueType, CancerType, ExhaustionGroup) 

cell_numbers <- data@meta.data %>%
    dplyr::count(pb_id, name = "n_cells")

pb_meta <- left_join(pb_meta, cell_numbers, by = "pb_id")

head(colnames(pb_counts), 10)
head(pb_meta$pb_id, 10)

pb_meta$pb_id <- as.character(pb_meta$pb_id)
pb_meta$pb_id <- gsub("_", "-", pb_meta$pb_id)

head(colnames(pb_counts), 10)
head(pb_meta$pb_id, 10)

pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$pb_id), ]
all(pb_meta$pb_id == colnames(pb_counts))

# Keep samples with more than 20 cells 
keep <- pb_meta$n_cells >= 20
table(keep)
pb_meta <- pb_meta[keep, ]
pb_counts <- pb_counts[, keep]

# Keep samples with more than 100 counts
keep_samples <- colSums(pb_counts) > 100
table(keep_samples)
pb_counts <- pb_counts[, keep_samples, drop = FALSE]
pb_meta <- pb_meta[keep_samples, , drop = FALSE]

rownames(pb_meta) <- pb_meta$pb_id

identical(colnames(pb_counts), rownames(pb_meta))

dge <- DGEList(counts = pb_counts)


# Compare between exhaustion groupings
table(pb_meta$ExhaustionGroup)
pb_meta$ExhaustionGroup <- factor(pb_meta$ExhaustionGroup, levels = c('Non-Exhausted','Precursor Exhausted', 'Exhausted'))

design <- model.matrix(~ 0 + ExhaustionGroup, data = pb_meta)

colnames(design) <- make.names(levels(pb_meta$ExhaustionGroup))
colnames(design)
design


# Filter low-expression genes
keep <- filterByExpr(dge, design)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalization
dge <- calcNormFactors(dge)

# voom
v <- voom(dge, design, plot = TRUE)

# Fit all three groups
fit <- lmFit(v, design)
colnames(design)
contrasts <- makeContrasts(
    PreEx = Precursor.Exhausted - Non.Exhausted,
    Ex = Exhausted - Non.Exhausted,
    PreEx_Ex       = Exhausted - Precursor.Exhausted,
    levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)



gene <- "HAVCR2"

# Get logFC and padj 

preex_res <- topTable(
    fit2, coef = "PreEx",
    number = Inf,
    sort.by = "none"
)[gene, ]

ex_res <- topTable(
    fit2, coef = "Ex",
    number = Inf,
    sort.by = "none"
)[gene, ]

preex_ex_res <- topTable(
    fit2, coef = "PreEx_Ex",
    number = Inf,
    sort.by = "none"
)[gene, ]

stats <- data.frame(
    group1 = c(
        "Non-Exhausted",
        "Non-Exhausted",
        "Precursor Exhausted"
    ),
    group2 = c(
        "Precursor Exhausted",
        "Exhausted",
        "Exhausted"
    ),
    logFC = c(
        preex_res$logFC,
        ex_res$logFC,
        preex_ex_res$logFC
    ),
    padj = c(
        preex_res$adj.P.Val,
        ex_res$adj.P.Val,
        preex_ex_res$adj.P.Val
    ),
    p = c(
      preex_res$P.Value,
      ex_res$P.Value,
      preex_ex_res$P.Value
    )
)

stats$p.adj.BH <- p.adjust(stats$p, method = 'BH')
stats


plot_df <- data.frame(
    sample = colnames(v$E),
    group = pb_meta$ExhaustionGroup,
    expression = as.numeric(v$E[gene, ])
)

plot_df$group <- factor(plot_df$group, levels = c("Non-Exhausted", "Precursor Exhausted","Exhausted"))

stats$label <- cut(
    stats$p.adj.BH,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
)
stats$number <- format(round(as.numeric(stats$p.adj.BH), 2), nsmall = 2)
stats$finallabel <- ifelse(stats$label == 'ns', stats$number, as.character(stats$label))

stats

ymax <- max(plot_df$expression, na.rm = TRUE)
yrange <- diff(range(plot_df$expression, na.rm = TRUE))

stats$y.position <- ymax + c(
    0.10,
    0.2,
    0.3
) * yrange

ggplot(plot_df,aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA,width = 0.55, alpha = 0.6) +
    geom_jitter(width = 0.10,size = 1,alpha = 0.8) +
    stat_pvalue_manual(
        stats,
        label = "finallabel", size = 4,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4
    ) +
    labs(x = NULL,y = "voom log2 HAVCR2 expression") +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none', 
          axis.text.x = element_text(colour = 'black', size = 10),
          axis.title.y = element_text(size = 14))
ggsave(paste0(output_path, 'TIM3_PB_DE_limma_exhaustion.pdf'), height = 5, width = 5)


# Q2. Is HAVCR2 expression different between primary and metastatic tumours? ----

data$pb_id <- interaction(
    data$batch,
    data$TissueType,
    data$CancerType,
    drop = T,
    sep = '--'
)
table(data$pb_id)
length(unique(data$pb_id))

pb_counts <- AggregateExpression(
    data,
    assays = 'RNA',
    slot = 'counts',
    group.by = 'pb_id',
    return.seurat = F
)$RNA %>%
    round()


pb_meta <- data@meta.data %>%
    as.data.frame() %>%
    distinct(pb_id, batch, TissueType, CancerType, ExhaustionGroup) 

cell_numbers <- data@meta.data %>%
    dplyr::count(pb_id, name = "n_cells")

pb_meta <- left_join(pb_meta, cell_numbers, by = "pb_id")

head(colnames(pb_counts), 10)
head(pb_meta$pb_id, 10)

pb_meta$pb_id <- as.character(pb_meta$pb_id)
pb_meta$pb_id <- gsub("_", "-", pb_meta$pb_id)

head(colnames(pb_counts), 10)
head(pb_meta$pb_id, 10)

pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$pb_id), ]
all(pb_meta$pb_id == colnames(pb_counts))

# Keep samples with more than 30 cells 
keep <- pb_meta$n_cells >= 30
table(keep)
pb_meta <- pb_meta[keep, ]
pb_counts <- pb_counts[, keep]


rownames(pb_meta) <- pb_meta$pb_id

identical(colnames(pb_counts), rownames(pb_meta))

dge <- DGEList(counts = pb_counts)


# Compare between exhaustion groupings
table(pb_meta$TissueType)
pb_meta$TissueType <- factor(pb_meta$TissueType, levels = c('Healthy Donor','Adjacent Tissue', 'Primary Tumour','Metastatic Tumour'))

design <- model.matrix(~ 0 + TissueType, data = pb_meta)

colnames(design) <- make.names(levels(pb_meta$TissueType))
colnames(design)
design


# Filter low-expression genes
keep <- filterByExpr(dge, design)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalization
dge <- calcNormFactors(dge)

# voom
v <- voom(dge, design, plot = TRUE)

# Fit all three groups
fit <- lmFit(v, design)
colnames(design)
contrasts <- makeContrasts(
    Meta_HD = Metastatic.Tumour - Healthy.Donor,
    Adj_HD = Adjacent.Tissue - Healthy.Donor,
    Primary_HD = Primary.Tumour - Healthy.Donor,
    Meta_Adj       = Metastatic.Tumour - Adjacent.Tissue,
    Primary_Adj = Primary.Tumour - Adjacent.Tissue,
    Meta_Primary = Metastatic.Tumour - Primary.Tumour,
    levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)



gene <- "HAVCR2"

# Get logFC and padj 

metahd_res <- topTable(
    fit2, coef = "Meta_HD",
    number = Inf,
    sort.by = "none"
)[gene, ]

adjhd_res <- topTable(
  fit2, coef = "Adj_HD",
  number = Inf,
  sort.by = "none"
)[gene, ]

primaryhd_res <- topTable(
    fit2, coef = "Primary_HD",
    number = Inf,
    sort.by = "none"
)[gene, ]

metaadj_res <- topTable(
    fit2, coef = "Meta_Adj",
    number = Inf,
    sort.by = "none"
)[gene, ]

primaryadj_res <- topTable(
    fit2, coef = "Primary_Adj",
    number = Inf,
    sort.by = "none"
)[gene, ]

metaprimary_res <- topTable(
    fit2, coef = "Meta_Primary",
    number = Inf,
    sort.by = "none"
)[gene, ]

stats <- data.frame(
    group1 = c(
        "Healthy Donor",
        "Healthy Donor",
        "Healthy Donor",
        "Adjacent Tissue",
        "Adjacent Tissue",
        "Primary Tumour"
    ),
    group2 = c(
      "Adjacent Tissue",
      "Primary Tumour",
      "Metastatic Tumour",
      "Primary Tumour",
      "Metastatic Tumour",
      "Metastatic Tumour"
    ),
    logFC = c(
        adjhd_res$logFC,
        primaryhd_res$logFC,
        metahd_res$logFC,
        primaryadj_res$logFC,
        metaadj_res$logFC,
        metaprimary_res$logFC
    ),
    padj = c(
      adjhd_res$adj.P.Val,
        primaryhd_res$adj.P.Val,
        metahd_res$adj.P.Val,
        primaryadj_res$adj.P.Val,
        metaadj_res$adj.P.Val,
        metaprimary_res$adj.P.Val
    ),
    p = c(
      adjhd_res$P.Value,
      primaryhd_res$P.Value,
      metahd_res$P.Value,
      primaryadj_res$P.Value,
      metaadj_res$P.Value,
      metaprimary_res$P.Value
    )
)

stats$p.adj.BH <- p.adjust(stats$p, method = 'BH')
stats



stats$number <- format(round(as.numeric(stats$p.adj.BH), 2), nsmall = 2)

stats$label <- cut(
    stats$p.adj.BH,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
)

stats$finallabel <- ifelse(stats$label == 'ns', stats$number, as.character(stats$label))

stats

plot_df <- data.frame(
    sample = colnames(v$E),
    group = pb_meta$TissueType,
    expression = as.numeric(v$E[gene, ])
)

plot_df$group <- factor(plot_df$group, levels = c('Healthy Donor','Adjacent Tissue', 'Primary Tumour','Metastatic Tumour'))


ymax <- max(plot_df$expression, na.rm = TRUE)
yrange <- diff(range(plot_df$expression, na.rm = TRUE))

stats$y.position <- ymax + c(
    0.10,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6
) * yrange

ggplot(plot_df,aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA,width = 0.55, alpha = 0.6) +
    geom_jitter(width = 0.10,size = 1,alpha = 0.8) +
    stat_pvalue_manual(
        stats,
        label = 'finallabel', size = 4,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4
    ) +
    labs(x = NULL,y = "voom log2 HAVCR2 expression") +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none', 
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(size = 14))
ggsave(paste0(output_path, 'TIM3_PB_DE_limma_tissue.pdf'), height = 5, width = 5)




# Q3. Is %HAVCR2+ different in any subsets or tissue types? ----


data$HAVCR2_pos <- ifelse(data@assays$RNA@counts['HAVCR2',]>0, 'Yes','No')
table(data$HAVCR2_pos)
table(data$HAVCR2_pos, data$TissueType)
table(data$HAVCR2_pos, data$simplecluster)

metadata <- data@meta.data %>%
    as.data.frame() %>%
    select(ExhaustionGroup, HAVCR2_pos, TissueType, batch, CancerType) %>%
    mutate(group = str_c(batch, CancerType, sep = '.')) %>%
    select(ExhaustionGroup, HAVCR2_pos, group, TissueType)

HAVCR2_pos_tissue <- metadata %>%
    group_by(group, TissueType) %>%
    summarise(pct.TIM3 = sum(HAVCR2_pos == 'Yes')*100/n(),
              n=n(),
              total.TIM3 = sum(HAVCR2_pos == 'Yes'),
              .groups = 'drop'
    ) %>%
    arrange(group, TissueType)


HAVCR2_pos_tissue$TissueType <- factor(HAVCR2_pos_tissue$TissueType , levels = c('Healthy Donor','Adjacent Tissue','Primary Tumour','Metastatic Tumour'))

# Non-paired test (wilcoxon rank-sum)
stat.test <- HAVCR2_pos_tissue %>%
    pairwise_wilcox_test(pct.TIM3 ~ TissueType) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "TissueType") 

stat.test$final_label = ifelse(stat.test$p.adj.signif == 'ns', format(round(as.numeric(stat.test$p.adj), 2), nsmall = 2), as.character(stat.test$p.adj.signif))

stat.test

ggplot(HAVCR2_pos_tissue, aes(x=TissueType, y = pct.TIM3, fill = TissueType)) +
    geom_boxplot(outlier.shape = NA,width = 0.55, alpha = 0.6) +
    geom_jitter(width = 0.10,size = 1,alpha = 0.8) +
    labs(x = NULL,y = "%HAVCR2+") +
    stat_pvalue_manual(
        stat.test,
        label = "final_label",
        tip.length = 0.01
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none', 
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(size = 14))
ggsave(paste0(output_path, 'TIM3_percent_tissue.pdf'), height = 5, width = 5)




# By exhaustion
HAVCR2_pos_exhaustion <- metadata %>%
    group_by(group, ExhaustionGroup) %>%
    summarise(pct.TIM3 = sum(HAVCR2_pos == 'Yes')*100/n(),
              n=n(),
              total.TIM3 = sum(HAVCR2_pos == 'Yes'),
              .groups = 'drop'
    ) %>%
    arrange(group, ExhaustionGroup)

HAVCR2_pos_exhaustion$ExhaustionGroup <- factor(HAVCR2_pos_exhaustion$ExhaustionGroup , levels = c('Non-Exhausted','Precursor Exhausted','Exhausted'))


df <- HAVCR2_pos_exhaustion %>%
    select(group, ExhaustionGroup, pct.TIM3)

pairs <- combn(levels(df$ExhaustionGroup), 2, simplify = FALSE)

# Paired tests between donors (wilcoxon signed-rank)
pairwise_tests <- map_df(pairs, function(pr) {
    tmp <- df %>%
        filter(ExhaustionGroup %in% pr) %>%
        pivot_wider(names_from = ExhaustionGroup, values_from = pct.TIM3) %>%
        drop_na(all_of(pr))
    
    wt <- wilcox.test(
        tmp[[pr[1]]],
        tmp[[pr[2]]],
        paired = TRUE,
        exact = FALSE
    )
    
    tibble(
        group1 = pr[1],
        group2 = pr[2],
        p.value = wt$p.value
    )
})



ymax <- max(df$pct.TIM3, na.rm = TRUE)

pairwise_tests <- pairwise_tests %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        y.position = ymax + seq_len(n()) * (0.08 * ymax + 1),
        label = cut(
            p.adj,
            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
            labels = c("***", "**", "*", "ns")
        ),
        final_label = ifelse(label == 'ns', format(round(as.numeric(p.adj), 2), nsmall = 2), as.character(label))
    )




ggplot(HAVCR2_pos_exhaustion, aes(x=ExhaustionGroup, y = pct.TIM3, fill = ExhaustionGroup)) +
    geom_boxplot(outlier.shape = NA,width = 0.55, alpha = 0.6) +
    geom_jitter(width = 0.10,size = 1,alpha = 0.8) +
    labs(x = NULL,y = "%HAVCR2+") +
    stat_pvalue_manual(
        pairwise_tests,
        label = "final_label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4,
        inherit.aes = FALSE) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none', 
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(size = 14))
ggsave(paste0(output_path, 'TIM3_percent_exhaustion.pdf'), height = 5, width = 5)

  




# Functional correlation analysis ----

functions_path <- paste0(output_path, 'functional_correlation/')
dir.create(functions_path, showWarnings = F)

# Extract basic data for analysis
gene_expression <- GetAssayData(data, assay = "RNA", layer = "data")[gene, ]
cell_types <- data@meta.data[[cell_cluster_name]]

# Get UMAP coordinates
umap_coords <- data@reductions$umap@cell.embeddings

# Create basic analysis dataframe
analysis_df <- data.frame(
    cell_id = colnames(data),
    cell_type = cell_types,
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    gene_expression = gene_expression,
    stringsAsFactors = FALSE
)

all_signatures <- grep('_Score1',colnames(data@meta.data), value = T)

# Extract scores and add to analysis dataframe
for (i in 1:length(all_signatures)) {
    signature_name <- all_signatures[i]
    if (signature_name %in% colnames(data@meta.data)) {
        analysis_df[[signature_name]] <- data@meta.data[[signature_name]]
    }
}

all_signatures <- gsub(x= all_signatures, pattern = '_Score1', replacement = '')

colnames(analysis_df) <- gsub(x= colnames(analysis_df), pattern = '_Score1', replacement = '')

# Function to calculate correlations
calculate_correlations <- function(signatures) {
    correlations_df <- data.frame()
    
    for (i in 1:length(signatures)) {
        score_col <- signatures[i]
        if (score_col %in% colnames(analysis_df)) {
            # Calculate correlation
            correlation <- cor(analysis_df$gene_expression, analysis_df[[score_col]], 
                               use = "complete.obs", method = "spearman")
            
            # Calculate p-value
            cor_test <- cor.test(analysis_df$gene_expression, analysis_df[[score_col]])
            p_value <- cor_test$p.value
            
            correlations_df <- rbind(correlations_df,
                                     data.frame(
                                         Signature = score_col,
                                         Correlation = correlation,
                                         P_value = p_value,
                                         stringsAsFactors = FALSE
                                     ))
        }
    }
    
    return(correlations_df)
}

# Calculate correlations for effector functions
all_correlations <- calculate_correlations(all_signatures)


# Adjust p-values for multiple testing
all_correlations$FDR <- p.adjust(all_correlations$P_value, method = "fdr")
all_correlations$Significant <- all_correlations$FDR < 0.01

# Sort by correlation strength
all_correlations <- all_correlations[order(abs(all_correlations$Correlation), decreasing = TRUE), ]

print("Top correlations:")
print(head(all_correlations, 10))


ggplot(all_correlations, aes(reorder(Signature, Correlation), y = Correlation, fill = Correlation)) +
    geom_col() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    coord_flip() +
    labs(x = NULL, y = 'Spearman Correlation Coefficient') +
    geom_hline(yintercept = 0, colour = "black") +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = 'black'))
ggsave(paste0(functions_path, gene, "_signaturecorrelations.pdf"), width = 5, height = 5)


print("Creating scatter plots for top correlations...")

all_correlations$Signature
signature_name <- 'Cytotoxicity'


ggplot(analysis_df, aes(x = gene_expression, y = get(signature_name))) +
    geom_point(alpha = 0.6, size = 0.1, color = 'grey30') +
    geom_smooth(method = "lm", color = "red", size = 0.8, se = F) +
    theme_bw() +
    labs(x = paste0(gene, " Expression"), y = paste(signature_name, "Score")) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme(legend.position = "bottom")

ggsave(paste0(functions_path, gene, "_vs_", tolower(signature_name), "_scatter.pdf"), width = 6, height = 5)









