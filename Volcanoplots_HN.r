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
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

# Get sample info ---------------------------------------------------------

samples <- c("HN200519A", "HN021219A", "HN120819A", "HN230620A")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

#Create figure directories:
if(!dir.exists(str_glue('{fig_dir}/{sample_name}'))){
  dir.create(str_glue('{fig_dir}/{sample_name}'))
}


#Load data for all ex and nonex
res<-read.csv(str_glue('{data_dir}/{sample_name}/{sample_name}_ex_v_non_ex.csv', header = TRUE))


#Basic plot, p value is 10e-6
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_basic_volcano.pdf"), width=11.6, height=8.2)
EnhancedVolcano(res, lab = res$X, x = "avg_log2FC", y = "p_val_adj") 
dev.off()


#Modify cut-offs for log2FC and P value; specify title; adjust point and label size - here p value is 10e-32
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_modified_volcano.pdf"), width=11.6, height=8.2)
 EnhancedVolcano(res,
    lab = res$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-6,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 6.0)
dev.off()

#Adjust colour and alpha for point shading p value is 10e-16
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_adjusted_colour_volcano.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res,
    lab = res$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-6,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 6.0,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1)
dev.off()

#Adjust p-value cut-off - here p value is 10e-5
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_all_adjusted_lines_volcano.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res,
    lab = res$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    xlim = c(-6, 6),
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
    gridlines.minor = FALSE)
dev.off()

#Load data for sentinel only
res<-read.csv(str_glue('{data_dir}/{sample_name}/{sample_name}_ex_sentinel.csv', header = TRUE))


#Basic plot, p value is 10e-6
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_sentinel_basic_volcano.pdf"), width=11.6, height=8.2)
EnhancedVolcano(res, lab = res$X, x = "avg_log2FC", y = "p_val_adj") 
dev.off()


#Modify cut-offs for log2FC and P value; specify title; adjust point and label size - here p value is 10e-32
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_sentinel_modified_volcano.pdf"), width=11.6, height=8.2)
 EnhancedVolcano(res,
    lab = res$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-6,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 6.0)
dev.off()

#Adjust colour and alpha for point shading p value is 10e-16
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_sentinel_adjusted_colour_volcano.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res,
    lab = res$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Expanding v Non-expanding clones',
    pCutoff = 10e-6,
    FCcutoff = 1.5,
    pointSize = 3.0,
    labSize = 6.0,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1)
dev.off()

#Adjust p-value cut-off - here p value is 10e-5
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_sentinel_adjusted_lines_volcano.pdf"), width=11.6, height=8.2)
  EnhancedVolcano(res,
    lab = res$X,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    xlim = c(-6, 6),
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
    gridlines.minor = FALSE)
dev.off()