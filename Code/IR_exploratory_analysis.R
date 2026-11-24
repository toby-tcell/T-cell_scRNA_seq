# Gene group analysis


# Libraries ----

library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ComplexHeatmap)


# Variables ----


input_path <- 'Data/CD8_processed.rds'
plot_path <- 'IR_plots/'
cell_cluster_col <- 'simplecluster'
dir.create(plot_path, showWarnings = F)

data <- readRDS(input_path)
gene_group <- c('PDCD1','LAG3','HAVCR2','TIGIT','CTLA4','BTLA')


# Dot Plots

DotPlot(data, features = gene_group, group.by = cell_cluster_col) +
    labs(y=NULL, x = NULL)+
    scale_x_discrete(guide = guide_axis(angle = 45)) 
ggsave(paste0(plot_path, 'dotplot.pdf'), height = 6, width = 9)


DotPlot(data, features = gene_group, group.by = 'TissueType') +
    labs(y=NULL, x = NULL)+
    scale_x_discrete(guide = guide_axis(angle = 45)) 
ggsave(paste0(plot_path, 'dotplot_tissue.pdf'), height = 6, width = 7)




