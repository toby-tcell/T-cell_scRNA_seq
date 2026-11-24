library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(ggrepel)
library(slingshot)
library(ggpubr)
library(colorspace)
library(ggnewscale)
library(RColorBrewer)
library(viridis)


# Loading correct data - ensure that data is CD8s only and annotated, for this pipeline


input_path <- 'Data/CD8_processed.rds'
cell_cluster_name <- 'simplecluster'
secondary_cluster_name <- 'TissueType'
pseudotime_plotpath <- 'Slingshot_Plots/'
data <- readRDS(input_path)

dir.create(pseudotime_plotpath, showWarnings = F)

# Slingshot ----

sce <- as.SingleCellExperiment(data)

sce_sling <- slingshot(sce,
                       clusterLabels = cell_cluster_name,
                       reducedDim = 'UMAP',
                       start.clus = "Naive",
                       shrink = T,
                       extend = 'n',
                       reweight = T,
                       approx_points = 150,
                       stretch = F)

saveRDS(sce_sling, "Data/scesling.rds")
sce_sling <- readRDS('Data/scesling.rds')

# Extract trajectory information
curves <- slingCurves(sce_sling)
lineages <- slingLineages(sce_sling)
pseudotime_matrix <- slingPseudotime(sce_sling)
curve_weights <- slingCurveWeights(sce_sling)


n_lineages <- ncol(pseudotime_matrix)
print(paste("Slingshot identified", n_lineages, "lineages"))
print("Lineages:")
for (i in 1:length(lineages)) {
    print(paste("Lineage", i, ":", paste(lineages[[i]], collapse = " -> ")))
}

print("Processing pseudotime and lineage assignments...")



umap_coords <- reducedDim(sce_sling, "UMAP")
df_plot <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    cluster = sce_sling$simplecluster,
    group = interaction(
        sce_sling$batch,
        sce_sling$CancerType
    )
)

curves <- slingCurves(sce_sling, as.df = TRUE) %>%
    arrange(Lineage, Order)



# For each cell, determine primary lineage (highest weight)
primary_lineage <- apply(curve_weights, 1, function(x) {
    if (all(is.na(x))) return(NA)
    which.max(x)
})

# Get primary pseudotime for each cell
primary_pseudotime <- sapply(1:ncol(sce_sling), function(i) {
    if (is.na(primary_lineage[i])) return(NA)
    pseudotime_matrix[i, primary_lineage[i]]
})

# Create lineage assignment
lineage_assignment <- factor(primary_lineage, levels = 1:n_lineages)

print(paste("Primary lineage distribution:"))
print(table(lineage_assignment, useNA = "ifany"))
print("Creating analysis dataframe...")

# Get coordinates for the specified reduction method
coords <- reducedDim(sce_sling, 'UMAP')

# Create comprehensive analysis dataframe
df_plot <- df_plot %>%
    mutate(
    primary_pseudotime = primary_pseudotime,
    primary_lineage = lineage_assignment,
    HAVCR2 = logcounts(sce_sling)["HAVCR2", ],
    stringsAsFactors = FALSE
)

endpoints <- curves %>% 
    dplyr::group_by(Lineage) %>% 
    slice_tail(n = 1) %>%
    mutate(nudge_y = case_when(
        Lineage %in% c(4,3) ~ 0.3,
        Lineage %in% c(1,2) ~ -0.3,
        TRUE ~ 0
    )) %>%
    mutate(nudge_x = case_when(
        Lineage %in% c(3) ~ -0.3,
        Lineage %in% c(2) ~ 0.3,
        TRUE ~ 0
    ))


ggplot(df_plot, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 0.5, alpha = 1) +
    geom_path(data = curves,
              aes(x = UMAP_1, y = UMAP_2, group = Lineage),
              colour = "black", linewidth = 1.2,
              arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
    geom_label(data = endpoints,
               aes(label = Lineage, x=UMAP_1, y = UMAP_2), colour = 'black', fill = 'white', fontface = 'bold',
               nudge_y = endpoints$nudge_y,
               nudge_x = endpoints$nudge_x) +
    theme(legend.position = 'none') +
    theme_classic() +
    labs(colour = "Cluster") 
    

ggsave(file.path(pseudotime_plotpath, 'UMAP_clusters.pdf'), height = 5, width = 7)


ggplot(df_plot, aes(x = UMAP1, y = UMAP2, colour = primary_pseudotime)) +
    geom_point(size = 0.5, alpha = 1) +
    geom_path(data = curves,
              aes(x = UMAP_1, y = UMAP_2, group = Lineage),
              colour = "black", linewidth = 1.2,
              arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
    geom_label(data = endpoints,
               aes(label = Lineage, x=UMAP_1, y = UMAP_2), colour = 'black', fill = 'white', fontface = 'bold',
               nudge_y = endpoints$nudge_y,
               nudge_x = endpoints$nudge_x) +
    theme(legend.position = 'none') +
    theme_classic() +
    scale_color_viridis() +
    labs(colour = "Pseudotime")

ggsave(file.path(pseudotime_plotpath, 'UMAP_pseudotime.pdf'), height = 5, width = 7)

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, colour = HAVCR2)) +
    geom_point(size = 0.5, alpha = 1) +
    geom_path(data = curves,
              aes(x = UMAP_1, y = UMAP_2, group = Lineage),
              colour = "black", linewidth = 1.2,
              arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
    geom_label(data = endpoints,
               aes(label = Lineage, x=UMAP_1, y = UMAP_2), colour = 'black', fill = 'white', fontface = 'bold',
               nudge_y = endpoints$nudge_y,
               nudge_x = endpoints$nudge_x) +
    theme(legend.position = 'none') +
    theme_classic() +
    scale_color_gradient(high = 'blue', low = 'grey') +
    labs(colour = "HAVCR2")

ggsave(file.path(pseudotime_plotpath, 'UMAP_TIM3.pdf'), height = 5, width = 7)


# TIM-3 vs pseudotime

# Overall TIM-3 vs primary pseudotime
valid_pseudotime <- !is.na(df_plot$primary_pseudotime)
ggplot(df_plot[valid_pseudotime, ], aes(x = primary_pseudotime, y = HAVCR2)) +
    geom_point(alpha = 0.6, size = 0.2, color = 'grey') +
    geom_smooth(method = 'gam', color = "red", size = 1, se = FALSE) +
    theme_bw() +
    labs(x = "Primary Pseudotime", y = paste("Normalized HAVCR2 Expression")) 

ggsave(paste0(pseudotime_plotpath, "TIM3_vs_primary_pseudotime.pdf"), width = 5, height = 4)


# TIM-3 vs pseudotime by lineage
valid_lineage <- !is.na(df_plot$primary_lineage) & !is.na(df_plot$primary_pseudotime)
ggplot(df_plot[valid_lineage, ], aes(x = primary_pseudotime, y = HAVCR2, color = primary_lineage)) +
    geom_point(alpha = 0.6, size = 0.2, color = 'grey') +
    geom_smooth(method = 'gam', size = 1, se = FALSE) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    theme_bw() +
    labs(x = "Primary Pseudotime", y = paste('Normalized HAVCR2 Expression'), color = "Primary Lineage") +
    theme(legend.position = "bottom")

ggsave(paste0(pseudotime_plotpath, "TIM3_vs_pseudotime_by_lineage.pdf"), width = 5, height = 4)


pt_quantiles <- c(0,5,10,15,20)
df_plot$pt_bin <- cut(df_plot$primary_pseudotime, breaks = pt_quantiles, include.lowest = TRUE)

table(df_plot$pt_bin)

bin_mapping <- c(
    "[0,5]"        = "0-5",
    "(5,10]"          = "5-10",
    "(10,15]"         = "10-15",
    "(15,20]"           = "15-20"
)

df_plot$bin <- unname(bin_mapping[as.character(df_plot$pt_bin)])
df_plot$bin <- factor(df_plot$bin, levels = c('0-5','5-10','10-15','15-20'))


# Plot pseudobulk TIM3 expression by pseudotime bins

df_pb <- df_plot %>%
    select(group, HAVCR2, bin) 

df_average <- df_pb %>%
    group_by(group, bin) %>%
    summarise(mean = mean(HAVCR2))

pairs <- list(
    c('0-5','15-20'),
    c('5-10','15-20'),
    c('10-15','15-20')
)

# Paired tests between donors
pairwise_tests <- map_df(pairs, function(pr) {
    tmp <- df_average %>%
        filter(bin %in% pr) %>%
        pivot_wider(names_from = bin, values_from = mean) %>%
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



ymax <- max(df_average$mean, na.rm = TRUE)

pairwise_tests <- pairwise_tests %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        y.position = ymax + seq_len(n()) * (0.001 * ymax + 0.2),
        label = cut(
            p.adj,
            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
            labels = c("***", "**", "*", "ns")
        ),
        final_label = ifelse(label == 'ns', format(round(as.numeric(p.adj), 2), nsmall = 2), as.character(label))
    )




ggplot(df_average, aes(x=bin, y = mean, fill = bin)) +
    geom_boxplot(outlier.shape = NA,width = 0.55, alpha = 0.6) +
    geom_jitter(width = 0.10,size = 1,alpha = 0.8) +
    labs(x = NULL,y = "Average TIM3 expression") +
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
    theme(legend.position = 'none') 
ggsave(paste0(pseudotime_plotpath, 'TIM3_average_bins.pdf'), height = 5, width = 5)

# Plot %TIM3 per group by pseudotime bins

df_pb <- df_plot %>%
    select(group, HAVCR2, bin) 

df_pct <- df_pb %>%
    mutate(HAVCR2_pos = ifelse(HAVCR2>0, 'Yes','No')) %>%
    group_by(group, bin) %>%
    summarise(pct.TIM3 = sum(HAVCR2_pos == 'Yes')*100/n(),
                  .groups = 'drop'
        ) 
    
pairs <- list(
    c('0-5','15-20'),
    c('5-10','15-20'),
    c('10-15','15-20')
)

# Paired tests between donors
pairwise_tests <- map_df(pairs, function(pr) {
    tmp <- df_pct %>%
        filter(bin %in% pr) %>%
        pivot_wider(names_from = bin, values_from = pct.TIM3) %>%
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



ymax <- max(df_pct$pct.TIM3, na.rm = TRUE)

pairwise_tests <- pairwise_tests %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        y.position = ymax + seq_len(n()) * (0.1 * ymax + 1),
        label = cut(
            p.adj,
            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
            labels = c("***", "**", "*", "ns")
        ),
        final_label = ifelse(label == 'ns', format(round(as.numeric(p.adj), 2), nsmall = 2), as.character(label))
    )




ggplot(df_pct, aes(x=bin, y = pct.TIM3, fill = bin)) +
    geom_boxplot(outlier.shape = NA,width = 0.55, alpha = 0.6) +
    geom_jitter(width = 0.10,size = 1,alpha = 0.8) +
    labs(x = NULL,y = "% TIM3+") +
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
    theme(legend.position = 'none') 
ggsave(paste0(pseudotime_plotpath, 'TIM3_percent_bins.pdf'), height = 5, width = 5)




# IR comparisons

df_irs <- df_plot %>%
    mutate(
        PDCD1 = logcounts(sce_sling)["PDCD1", ],
        LAG3 = logcounts(sce_sling)["LAG3", ],
        BTLA = logcounts(sce_sling)["BTLA", ],
        CTLA4 = logcounts(sce_sling)["CTLA4", ],
        TIGIT = logcounts(sce_sling)["TIGIT", ],
    ) %>%
    pivot_longer(cols = c((ncol(.)-4):ncol(.),'HAVCR2'), names_to = 'Gene', values_to = 'Expression') %>%
    select(Gene, Expression, primary_pseudotime, bin)

ggplot(df_irs, aes(x = primary_pseudotime, y = Expression, colour = Gene)) +
    geom_smooth(method = 'gam', size = 1, se = FALSE) +
    theme_bw() +
    labs(x = "Primary Pseudotime", y = paste("Normalized Expression"))

ggsave(paste0(pseudotime_plotpath, "IRs_vs_pseudotime.pdf"), width = 5, height = 3)


# Plot Gene expression by pseudotime bins



df_pct <- df_irs %>%
    group_by(Gene, bin) %>%
    summarise(n=sum(Expression)) %>%
    group_by(Gene) %>%
    mutate(pct.expression = prop.table(n)*100)


ggplot(df_pct, aes(x = Gene, y = pct.expression, fill = bin)) +
    geom_bar(stat = 'identity', show.legend = F, color = 'black', linewidth = 0.2) +
    ylab('%Expression') +
    geom_text(aes(label = ifelse(pct.expression < 3, '' , paste0(bin))), 
              position = position_stack(vjust = 0.5), size = 4) +
    labs(x = NULL, y = '%Expression in Bin') +
    scale_fill_discrete_diverging(palette = 'Blue-Red 2') +
    theme_bw()
ggsave(paste0(pseudotime_plotpath, 'Pct_bin_byIR.pdf'), height = 5, width = 5)



