#Do some example featureplots for the paper
#HN021219A, HN120819A, HN200519A, HN230620A

#The data is here:
#data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
#fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/featureplots"
#code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(Nebulosa)
library(glue)
library(SeuratDisk)
library(patchwork)
})

print("set directories")
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/featureplots"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

print("set sample names")
samples <- c("HN021219A", "HN120819A", "HN200519A", "HN230620A", "HN070219A", "HN170419A", "HN120520A")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("load in the data")
x<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_clones_final_immune_object.rds'))

print("change default assay")
DefaultAssay(x)<-"SCT"

print("do featureplot")
pdf(str_glue("{fig_dir}/{sample_name}_featureplot.pdf"), width=11.6, height=8.2)
FeaturePlot(x, features = c("EPCAM", "MS4A1", "CD4", "CD8B", "NCAM1", "CD14", "FCER1A", "SERPINF1"), ncol = 4, order = TRUE) & NoLegend() & NoAxes()
dev.off()
