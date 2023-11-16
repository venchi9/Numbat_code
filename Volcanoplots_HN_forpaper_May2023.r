#Try volcano plots for the HN data
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

# Load required R packages
suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(EnhancedVolcano)
})

#Define directories
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures/volcanoplots"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

# Get sample info ---------------------------------------------------------

samples <- c("HN200519A", "HN021219A", "HN120819A", "HN230620A", "combined")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

#Create figure directories:
if(!dir.exists(str_glue('{fig_dir}/volcanoplots/{sample_name}'))){
  dir.create(str_glue('{fig_dir}/volcanoplots/{sample_name}'))
}


#Load data for all ex and nonex
res<-read.csv(str_glue('{data_dir}/{sample_name}_DEG_list.csv', header = TRUE))

#load in new combined list from new object
res<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined_cancer_DEG.csv")
sample_name<-"combined"

#Filter for coding genes:
coding<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/master_coding_gene_list.csv")
coding_gene<-coding$gene

#Subset the res table for genes in "coding_gene list"
res_coding<-filter(res, res$X %in% coding_gene)

#Basic plot, p value is 10e-6
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_basic_volcano.pdf"), width=11.6, height=8.2)
EnhancedVolcano(res_coding, lab = res_coding$X, x = "avg_log2FC", y = "p_val_adj") 
dev.off()


#Modify cut-offs for log2FC and P value; specify title; adjust point and label size - here p value is 10e-32
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_modified_volcano_FC1_5.pdf"), width=11.6, height=8.2)
 EnhancedVolcano(res_coding,
    lab = res_coding$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-5,
    FCcutoff = 0,
    pointSize = 3.0,
    labSize = 6.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75)
dev.off()

#Adjust colour and alpha for point shading p value is 10e-16
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_adjusted_colour_volcano_FC1_5.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res_coding,
    lab = res_coding$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-4,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 6.0,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1,
    drawConnectors = TRUE,
    widthConnectors = 0.75)
dev.off()

#Adjust p-value cut-off - here p value is 10e-5
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_adjusted_lines_volcano_FC1_5.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res_coding,
    lab = res_coding$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    xlim = c(-4, 4),
    title = 'Expanding v Non-expanding clones',
    subtitle = paste0('p-value cutoff (red line) drawn ',
      'at equivalent of adjusted p=0.00001'),
    pCutoff = 0.00001,
    pCutoffCol = 'p_val_adj',
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 6.0,
    colAlpha = 1,
    cutoffLineType = 'solid',
    cutoffLineCol = 'red2',
    cutoffLineWidth = 2.5,
    hline = c(10e-20,
      10e-20 * 10e-30,
      10e-20 * 10e-60,
      10e-20 * 10e-90),
    hlineCol = c('black', 'black', 'black', 'black'),
    hlineType = c('longdash', 'longdash', 'dotdash', 'dotdash'),
    hlineWidth = c(0.4, 0.4, 0.8, 0.8),
    gridlines.major = FALSE,
    gridlines.minor = FALSE, 
    drawConnectors = TRUE,
    widthConnectors = 0.75)
dev.off()

#FC2


#Modify cut-offs for log2FC and P value; specify title; adjust point and label size - here p value is 10e-32
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_modified_volcano_FC2.pdf"), width=11.6, height=8.2)
 EnhancedVolcano(res_coding,
    lab = res_coding$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-4,
    FCcutoff = 2.0,
    pointSize = 3.0,
    labSize = 6.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75)
dev.off()

#Adjust colour and alpha for point shading p value is 10e-16
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_adjusted_colour_volcano_FC2.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res_coding,
    lab = res_coding$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-4,
    FCcutoff = 2,
    pointSize = 3.0,
    labSize = 6.0,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1,
    drawConnectors = TRUE,
    widthConnectors = 0.75)
dev.off()

#Adjust p-value cut-off - here p value is 10e-5
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_adjusted_lines_volcano_FC2.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res_coding,
    lab = res_coding$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    xlim = c(-4, 4),
    title = 'Expanding v Non-expanding clones',
    subtitle = paste0('p-value cutoff (red line) drawn ',
      'at equivalent of adjusted p=0.00001'),
    pCutoff = 0.00001,
    pCutoffCol = 'p_val_adj',
    FCcutoff = 2,
    pointSize = 3.0,
    labSize = 6.0,
    colAlpha = 1,
    cutoffLineType = 'solid',
    cutoffLineCol = 'red2',
    cutoffLineWidth = 2.5,
    hline = c(10e-20,
      10e-20 * 10e-30,
      10e-20 * 10e-60,
      10e-20 * 10e-90),
    hlineCol = c('black', 'black', 'black', 'black'),
    hlineType = c('longdash', 'longdash', 'dotdash', 'dotdash'),
    hlineWidth = c(0.4, 0.4, 0.8, 0.8),
    gridlines.major = FALSE,
    gridlines.minor = FALSE, 
    drawConnectors = TRUE,
    widthConnectors = 0.75)
dev.off()



