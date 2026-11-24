# TIM-3 Co-expression and Pathway Analysis for scRNA-seq T-cells
# Memory-efficient approach for large datasets (110K cells)

# ============================================================================
# SETUP AND LIBRARIES
# ============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(fgsea)
library(viridis)
library(cowplot)
library(ggrepel)


# ============================================================================
# DATA LOADING AND VALIDATION
# ============================================================================

# Load Seurat object - modify path as needed
seurat_obj <- readRDS("../Data/CD8s_Tex_manascore_final.rds")


# Validate data structure
cat("Seurat object summary:\n")
print(seurat_obj)
cat("\nMetadata columns:\n")
print(colnames(seurat_obj@meta.data))

# Check required columns
required_cols <- c("cell.type", "TissueType", "CancerType")
missing_cols <- required_cols[!required_cols %in% colnames(seurat_obj@meta.data)]
if (length(missing_cols) > 0) {
    stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
}

# Create derived condition columns
table(seurat_obj@meta.data$TissueType)
seurat_obj$condition <- ifelse(seurat_obj$TissueType %in% c("Healthy donor","Uninvolved normal tissue"), "Healthy", "Cancer")
seurat_obj$detailed_condition <- seurat_obj$TissueType

# Create patient_id if missing
table(seurat_obj@meta.data$batch)
colnames(seurat_obj)
seurat_obj$patient_id <- paste(seurat_obj$batch, seurat_obj$CancerType, sep = "_")

# Print data distribution
cat("\nData distribution:\n")
print(table(seurat_obj$condition))
print(table(seurat_obj$detailed_condition))
print(table(seurat_obj$patient_id))

# Create output directory
if (!dir.exists("plots")) dir.create("plots")

# ============================================================================
# TIM-3 EXPRESSION VALIDATION
# ============================================================================

tim3_gene <- "HAVCR2"
if (!tim3_gene %in% rownames(seurat_obj)) {
    stop("TIM-3 (HAVCR2) not found in dataset")
}

cat("\nTIM-3 expression summary:\n")
cat("Expressing cells:", sum(seurat_obj@assays$RNA@data[tim3_gene,] > 0), "/", ncol(seurat_obj), "\n")
cat("Mean expression:", round(mean(seurat_obj@assays$RNA@data[tim3_gene,]), 3), "\n")

# ============================================================================
# CORRELATION ANALYSIS FUNCTIONS
# ============================================================================

perform_correlation_analysis <- function(seurat_subset, condition_name) {
    cat("Analyzing", condition_name, "with", ncol(seurat_subset), "cells\n")
    
    # Get sparse expression matrix (DON'T convert to dense)
    expr_matrix <- seurat_subset@assays$RNA@data
    
    # Filter genes: expressed in at least 10% of cells
    min_cells <- ceiling(ncol(expr_matrix) * 0.1)
    gene_counts <- Matrix::rowSums(expr_matrix > 0)
    genes_keep <- names(gene_counts)[gene_counts >= min_cells]
    
    # Always keep TIM-3 even if it doesn't meet threshold
    if (!tim3_gene %in% genes_keep) {
        genes_keep <- c(genes_keep, tim3_gene)
        cat("TIM-3 added despite low expression in", condition_name, "\n")
    }
    expr_matrix <- expr_matrix[genes_keep, ]
    
    cat("Analyzing", length(genes_keep), "genes\n")
    
    # Extract TIM-3 expression as dense vector (only one gene)
    tim3_expr <- as.numeric(expr_matrix[tim3_gene, ])
    
    # Calculate correlations gene by gene (keeping sparse)
    correlations <- numeric(length(genes_keep))
    p_values <- numeric(length(genes_keep))
    names(correlations) <- genes_keep
    names(p_values) <- genes_keep
    
    # Use gene names instead of numeric indices
    for (gene in genes_keep) {
        gene_expr <- as.numeric(expr_matrix[gene, ])
        cor_result <- cor.test(tim3_expr, gene_expr, method = "spearman")
        correlations[gene] <- cor_result$estimate
        p_values[gene] <- cor_result$p.value
    }
    
    # Create results
    results <- data.frame(
        gene = names(correlations),
        correlation = correlations,
        p_value = p_values,
        condition = condition_name,
        fdr = p.adjust(p_values, method = "fdr"),
        stringsAsFactors = FALSE
    ) %>%
        arrange(desc(abs(correlation)))
    
    return(results)
}

calculate_differential_coexpression <- function(cancer_cor, healthy_cor) {
    common_genes <- intersect(cancer_cor$gene, healthy_cor$gene)
    
    diff_results <- data.frame(
        gene = common_genes,
        cancer_cor = cancer_cor$correlation[match(common_genes, cancer_cor$gene)],
        healthy_cor = healthy_cor$correlation[match(common_genes, healthy_cor$gene)]
    ) %>%
        mutate(
            cor_difference = cancer_cor - healthy_cor,
            cancer_z = 0.5 * log((1 + cancer_cor) / (1 - cancer_cor)),
            healthy_z = 0.5 * log((1 + healthy_cor) / (1 - healthy_cor))
        )
    
    # Statistical testing
    n_cancer <- sum(seurat_obj$condition == "Cancer")
    n_healthy <- sum(seurat_obj$condition == "Healthy")
    
    diff_results <- diff_results %>%
        mutate(
            se_diff = sqrt(1/(n_cancer-3) + 1/(n_healthy-3)),
            z_score = (cancer_z - healthy_z) / se_diff,
            p_value = 2 * pnorm(-abs(z_score)),
            fdr = p.adjust(p_value, method = "fdr")
        ) %>%
        arrange(fdr)
    
    return(diff_results)
}

calculate_cancerdifferential_coexpression <- function(cancer_cor, healthy_cor) {
    common_genes <- intersect(cancer_cor$gene, healthy_cor$gene)
    
    diff_results <- data.frame(
        gene = common_genes,
        cancer_cor = cancer_cor$correlation[match(common_genes, cancer_cor$gene)],
        healthy_cor = healthy_cor$correlation[match(common_genes, healthy_cor$gene)]
    ) %>%
        mutate(
            cor_difference = cancer_cor - healthy_cor,
            cancer_z = 0.5 * log((1 + cancer_cor) / (1 - cancer_cor)),
            healthy_z = 0.5 * log((1 + healthy_cor) / (1 - healthy_cor))
        )
    
    # Statistical testing
    n_cancer <- sum(seurat_obj$detailed_condition == "Primary tumor tissue")
    n_healthy <- sum(seurat_obj$detailed_condition == "Metastatic tumor tissue")
    
    diff_results <- diff_results %>%
        mutate(
            se_diff = sqrt(1/(n_cancer-3) + 1/(n_healthy-3)),
            z_score = (cancer_z - healthy_z) / se_diff,
            p_value = 2 * pnorm(-abs(z_score)),
            fdr = p.adjust(p_value, method = "fdr")
        ) %>%
        arrange(fdr)
    
    return(diff_results)
}
# ============================================================================
# MAIN CORRELATION ANALYSIS
# ============================================================================

cat("\n=== CORRELATION ANALYSIS ===\n")

# Primary analysis: Cancer vs Healthy
cancer_cells <- seurat_obj[, seurat_obj$condition == "Cancer"]
healthy_cells <- seurat_obj[, seurat_obj$condition == "Healthy"]

cancer_correlations <- perform_correlation_analysis(cancer_cells, "Cancer")
healthy_correlations <- perform_correlation_analysis(healthy_cells, "Healthy")

# Differential co-expression

diff_coexpr <- calculate_differential_coexpression(cancer_correlations, healthy_correlations)

diff_coexpr <- read.csv('Coexpression+pathways/tim3_differential_coexpression.csv')

# Tissue-specific analysis (if sufficient cells)
tissue_correlations <- list()
tissue_types <- unique(seurat_obj@meta.data$TissueType)

for (tissue in tissue_types) {
    cells <- seurat_obj[, seurat_obj$TissueType == tissue]
    if (ncol(cells) >= 100) {
        tissue_name <- gsub(" ", "_", gsub(" tissue", "", tissue))
        tissue_correlations[[tissue_name]] <- perform_correlation_analysis(cells, tissue_name)
    }
}
primary_correlations <- tissue_correlations[["Primary_tumor"]]
metastatic_correlations <- tissue_correlations[["Metastatic_tumor"]]


cancer_diff_coexpr <- calculate_cancerdifferential_coexpression(metastatic_correlations, primary_correlations)
# Clean up memory
rm(cancer_cells, healthy_cells)
gc()

# ============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# ============================================================================

cat("\n=== PATHWAY ANALYSIS ===\n")

# Load gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
immune_sets <- msigdbr(species = "Homo sapiens", category = "C7")

# Convert to list format for fgsea
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
immune_list <- split(immune_sets$gene_symbol, immune_sets$gs_name)

# Filter immune sets for T-cell specific pathways
tcell_patterns <- c("CD8", "TCELL", 
                    "EXHAUSTED", "NAIVE", "MEMORY", "EFFECTOR", "CYTOTOXIC")
tcell_names <- immune_sets$gs_name[grepl(paste(tcell_patterns, collapse = "|"), 
                                         immune_sets$gs_name, ignore.case = TRUE)]
tcell_list <- immune_list[names(immune_list) %in% tcell_names]

cat("Found", length(tcell_list), "T-cell specific pathways\n")

# Create custom T-cell signatures based on literature
custom_tcell_sets <- list(
    "TIM3_EXHAUSTION_CORE" = c("HAVCR2", "PDCD1", "LAG3", "TIGIT", "TOX", "CTLA4", 
                               "TOX2", "CD244", "BTLA", "CD160", "KLRG1", "CX3CR1"),
    "CYTOTOXIC_FUNCTION" = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", 
                             "IFNG", "TNF", "FASLG", "TRAIL", "CD107A", "NKG7"),
    "T_CELL_ACTIVATION" = c("CD69", "CD25", "IL2RA", "ICOS", "CD40LG", "CD44", 
                            "CD38", "HLA-DRA", "TNFRSF9", "TNFRSF4", "CD137", "OX40"),
    "T_CELL_MEMORY" = c("IL7R", "CCR7", "SELL", "TCF7", "BCL2", "CD127", 
                        "LDHB", "CPT1A", "FOXO1", "KLF2", "S1PR1", "CD28"),
    "METABOLIC_DYSFUNCTION" = c("LDHA", "HK2", "PFKP", "ENO1", "GAPDH", "PKM", 
                                "ALDOA", "TPI1", "PGAM1", "PGK1", "PKLR", "G6PD"),
    "TREG_SIGNATURE" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "IL10", 
                         "TIGIT", "ICOS", "HELIOS", "NRIP1", "LRRC32", "CD39")
)

# GSEA function
perform_gsea <- function(gene_rankings, gene_sets, set_name) {
    cat("Running GSEA for", set_name, "with", length(gene_sets), "pathways\n")
    
    results <- fgsea(
        pathways = gene_sets,
        stats = gene_rankings,
        minSize = 10,
        maxSize = 500
    )
    results$gene_set_collection <- set_name
    return(results)
}

# Prepare gene rankings (correlation values)
cancer_rankings <- setNames(cancer_correlations$correlation, cancer_correlations$gene)
healthy_rankings <- setNames(healthy_correlations$correlation, healthy_correlations$gene)

# Remove genes with NA correlations
cancer_rankings <- cancer_rankings[!is.na(cancer_rankings)]
healthy_rankings <- healthy_rankings[!is.na(healthy_rankings)]

# Run GSEA for all pathway collections
cat("Running comprehensive GSEA analysis...\n")

all_gsea_results <- rbind(
    # Hallmark pathways
    cbind(perform_gsea(cancer_rankings, hallmark_list, "Hallmark"), condition = "Cancer"),
    # cbind(perform_gsea(healthy_rankings, hallmark_list, "Hallmark"), condition = "Healthy"),
    
    # T-cell specific pathways from MSigDB
    cbind(perform_gsea(cancer_rankings, tcell_list, "T_Cell_MSigDB"), condition = "Cancer")
    # cbind(perform_gsea(healthy_rankings, tcell_list, "T_Cell_MSigDB"), condition = "Healthy"),
    
    # Custom T-cell signatures
    # cbind(perform_gsea(cancer_rankings, custom_tcell_sets, "T_Cell_Custom"), condition = "Cancer"),
    # cbind(perform_gsea(healthy_rankings, custom_tcell_sets, "T_Cell_Custom"), condition = "Healthy")
)

# Summary of significant results
cat("\nGSEA Results Summary:\n")
sig_results <- all_gsea_results %>%
    filter(padj < 0.05) %>%
    group_by(condition, gene_set_collection) %>%
    summarise(n_significant = n(), .groups = 'drop')

print(sig_results)

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

plot_tim3_expression <- function(save_plot = TRUE) {
    tim3_data <- FetchData(seurat_obj, vars = c(tim3_gene, "condition", "TissueType", "cell.type"))
    colnames(tim3_data) <- c("TIM3_expression", "condition", "tissue_type", "cell_type")
    
    p1 <- ggplot(tim3_data, aes(x = condition, y = TIM3_expression, fill = condition)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.2, fill = "grey", outlier.shape = NA, alpha = 0.5) +
        scale_fill_manual(values = c("Cancer" = "#d62728", "Healthy" = "#1f77b4")) +
        labs(title = "TIM-3 Expression: Cancer vs Healthy",
             x = "Condition", y = "TIM-3 Expression") +
        theme_bw() +
        theme(legend.position = "none")
    
    p2 <- ggplot(tim3_data, aes(x = tissue_type, y = TIM3_expression, fill = condition)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.2, fill = "grey", outlier.shape = NA, alpha = 0.5) +
        scale_fill_manual(values = c("Cancer" = "#d62728", "Healthy" = "#1f77b4")) +
        labs(title = "TIM-3 Expression by Tissue Type", x = "Tissue Type", y = "TIM-3 Expression") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    combined_plot <- plot_grid(p1, p2, ncol = 2)
    
    if (save_plot) {
        ggsave("plots/tim3_expression_comparison.pdf", combined_plot, width = 14, height = 6)
        cat("Saved: plots/tim3_expression_comparison.pdf\n")
    }
    return(combined_plot)
}

plot_correlations <- function(cor_results, condition_name, top_n = 30, save_plot = TRUE) {
    plot_data <- rbind(
        cor_results %>% filter(fdr < 0.05, correlation > 0) %>% head(top_n/2),
        cor_results %>% filter(fdr < 0.05, correlation < 0) %>% head(top_n/2)
    )
    
    p <- ggplot(plot_data, aes(x = reorder(gene, correlation), y = correlation)) +
        geom_col(aes(fill = correlation > 0)) +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4")) +
        labs(title = paste("Top TIM-3 Correlations -", condition_name),
             x = "Gene", y = "Spearman Correlation") +
        theme_bw() +
        theme(legend.position = "none")
    
    if (save_plot) {
        filename <- paste0("plots/correlations_", tolower(gsub(" ", "_", condition_name)), ".pdf")
        ggsave(filename, p, width = 10, height = 8)
        cat("Saved:", filename, "\n")
    }
    return(p)
}

plot_differential_coexpression <- function(diff_results) {
    plot_data <- diff_results %>%
        filter(fdr < 0.05) %>%
        arrange(desc(abs(cor_difference))) %>%
        mutate(fdr_adj = pmax(fdr, 1e-50),
               neg_log_fdr = -log10(fdr_adj))
    
    plot_data$label <- ifelse(1:nrow(plot_data) <= 20, plot_data$gene, "")
    
    p <- ggplot(plot_data, aes(x = healthy_cor, y = cancer_cor)) +
        geom_point(aes(color = cor_difference, size = neg_log_fdr), alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        scale_size_continuous(range = c(0, 10)) +
        geom_text_repel(aes(label = label), 
                        size = 4,
                        max.overlaps = 20,       
                        box.padding = 0.3,       
                        point.padding = 0.4,     
                        segment.size = 0.25,       
                        segment.alpha = 0.6,     
                        force = 1,                
                        force_pull = 3,          
                        min.segment.length = 0.1, 
                        seed = 42) +
        scale_color_gradient2(low = "#1f77b4", mid = "white", high = "#d62728",
                              name = "Correlation\nDifference") +
        labs(title = "Differential Co-expression with TIM-3", subtitle = 'with top 20 differentially co-expressed genes',
             x = "Healthy Correlation", y = "Cancer Correlation") +
        theme_bw()
    
    ggsave(paste0("Coexpression+Pathways/plots/differential_coexpression.pdf"), p, width = 15, height = 12)
    return(p)
}

plot_cancer_differential_coexpression <- function(diff_results, save_plot = TRUE) {
    plot_data <- diff_results %>%
        filter(fdr < 0.05) %>%
        arrange(desc(abs(cor_difference))) %>%
        mutate(fdr_adj = pmax(fdr, 1e-50),  
           neg_log_fdr = -log10(fdr_adj))
    
    plot_data$label <- ifelse(1:nrow(plot_data) <= 25, plot_data$gene, "")
    
    
    
    p <- ggplot(plot_data, aes(x = healthy_cor, y = cancer_cor)) +
        geom_point(aes(color = cor_difference, 
                       size = neg_log_fdr), 
                   alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        geom_text_repel(aes(label = label), 
                        size = 1.5,
                        max.overlaps = 20,       
                        box.padding = 0.35,       
                        point.padding = 0.5,     
                        segment.size = 0.2,       
                        segment.alpha = 0.6,     
                        force = 1,                
                        force_pull = 1,          
                        min.segment.length = 0.1, 
                        seed = 42) +
        scale_color_gradient2(low = "#1f77b4", mid = "white", high = "#d62728",
                              name = "Correlation\nDifference") +
        labs(title = "Differential Co-expression with TIM-3",
             x = "Primary Tumour Correlation", y = "Metastatic Tumour Correlation") +
        theme_bw()
    
    if (save_plot) {
        ggsave(paste0("plots/cancer_differential_coexpression.pdf"), p, width = 10, height = 8)
    }
    return(p)
}

plot_gsea_results <- function(gsea_results, pathway_collection = "Hallmark", save_plot = TRUE, max_pathways = 20) {
    # Get all pathways that are significant in either condition
    significant_pathways <- gsea_results %>%
        filter(gene_set_collection == pathway_collection, padj < 0.05) %>%
        pull(pathway) %>%
        unique()
    
    if (length(significant_pathways) == 0) {
        cat("No significant pathways for", pathway_collection, "\n")
        return(NULL)
    }
    
    # Get data for ALL conditions for these pathways (including non-significant)
    plot_data <- gsea_results %>%
        filter(gene_set_collection == pathway_collection, 
               pathway %in% significant_pathways) %>%
        mutate(pathway_clean = gsub("HALLMARK_|_", " ", pathway)) %>%
        mutate(pathway_clean = stringr::str_to_title(pathway_clean)) %>%
        # Mark non-significant results
        mutate(significant = padj < 0.05,
               NES_display = NES) %>%  # Keep NES even if not significant
        group_by(pathway) %>%
        mutate(max_abs_nes = max(abs(NES))) %>%  # Calculate max_abs_nes within the plot_data
        ungroup()
    
    # Limit to top pathways based on maximum absolute NES across conditions
    top_pathways <- plot_data %>%
        select(pathway, max_abs_nes) %>%
        distinct() %>%
        arrange(desc(max_abs_nes)) %>%
        head(max_pathways) %>%
        pull(pathway)
    
    plot_data <- plot_data %>%
        filter(pathway %in% top_pathways)
    
    # Shorten pathway names if too long
    plot_data$pathway_short <- ifelse(nchar(plot_data$pathway_clean) > 35, 
                                      paste0(substr(plot_data$pathway_clean, 1, 32), "..."), 
                                      plot_data$pathway_clean)
    
    p <- ggplot(plot_data, aes(x = reorder(pathway_short, max_abs_nes), y = NES_display, fill = condition)) +
        geom_col(position = "dodge", width = 0.7, 
                 aes(alpha = significant)) +  # Make non-significant bars more transparent
        coord_flip() +
        scale_fill_manual(values = c("Cancer" = "#d62728", "Healthy" = "#1f77b4")) +
        scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.4), guide = "none") +
        labs(title = paste("Top", pathway_collection, "Pathways"),
             subtitle = "Significant pathways (FDR < 0.05) shown in full color, non-significant in faded color",
             x = "Pathway", y = "Normalized Enrichment Score") +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 9),
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 10),
            legend.position = "bottom"
        )
    
    if (save_plot) {
        filename <- paste0("plots/gsea_", tolower(pathway_collection), ".pdf")
        plot_height <- max(6, min(12, length(top_pathways) * 0.4 + 3))
        ggsave(filename, p, width = 10, height = plot_height)
        cat("Saved:", filename, "\n")
    }
    return(p)
}


# ============================================================================
# GENERATE ALL VISUALIZATIONS
# ============================================================================

cat("\n=== GENERATING PLOTS ===\n")

plot_tim3_expression()
plot_correlations(cancer_correlations, "Cancer")
plot_correlations(healthy_correlations, "Healthy")

# Plot tissue-specific results
for (name in names(tissue_correlations)) {
    plot_correlations(tissue_correlations[[name]], name)
}

plot_differential_coexpression(diff_coexpr)
plot_cancer_differential_coexpression(cancer_diff_coexpr)

# Plot GSEA results
for (collection in unique(all_gsea_results$gene_set_collection)) {
    plot_gsea_results(all_gsea_results, collection)
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Main results
write.csv(cancer_correlations, "tim3_cancer_correlations.csv", row.names = FALSE)
write.csv(healthy_correlations, "tim3_healthy_correlations.csv", row.names = FALSE)
write.csv(diff_coexpr, "tim3_differential_coexpression.csv", row.names = FALSE)
write.csv(cancer_diff_coexpr, "tim3_cancer_differential_coexpression.csv", row.names = FALSE)


# Tissue-specific results
for (name in names(tissue_correlations)) {
    filename <- paste0("tim3_", tolower(name), "_correlations.csv")
    write.csv(tissue_correlations[[name]], filename, row.names = FALSE)
}

# Summary statistics
summary_stats <- data.frame(
    Analysis = c("Cancer correlations", "Healthy correlations", "Differential coexpression", 
                 "Cancer pathways", "Healthy pathways"),
    Significant_Results = c(
        sum(cancer_correlations$fdr < 0.05),
        sum(healthy_correlations$fdr < 0.05),
        sum(diff_coexpr$fdr < 0.05),
        sum(all_gsea_results$condition == "Cancer" & all_gsea_results$padj < 0.05),
        sum(all_gsea_results$condition == "Healthy" & all_gsea_results$padj < 0.05)
    )
)
write.csv(summary_stats, "analysis_summary.csv", row.names = FALSE)

# ============================================================================
# RESULTS SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Dataset:", ncol(seurat_obj), "T-cells\n")
cat("Cancer cells:", sum(seurat_obj$condition == "Cancer"), "\n")
cat("Healthy cells:", sum(seurat_obj$condition == "Healthy"), "\n")

cat("\nSignificant results (FDR < 0.05):\n")
print(summary_stats)

cat("\nTop 5 genes positively correlated with TIM-3 in cancer:\n")
top_cancer <- cancer_correlations %>% 
    filter(fdr < 0.05, correlation > 0) %>% 
    head(5) %>%
    select(gene, correlation, fdr)
print(top_cancer)

cat("\nTop 5 differential co-expression genes:\n")
top_diff <- diff_coexpr %>%
    filter(fdr < 0.05) %>%
    head(5) %>%
    select(gene, cancer_cor, healthy_cor, cor_difference, fdr)
print(top_diff)

cat("\nAll results saved to CSV files and plots to 'plots/' directory\n")

# Save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")