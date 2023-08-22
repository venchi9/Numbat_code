#Need to take all the pooled objects normalise them and generate plots so I can roughly annotate them.

# Load required R packages
suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(glue)
})

#Define directories
data_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures"
code_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"


samples <- c("HN021219A", "HN070219A", "HN120520A","HN120819A", "HN200519A", "HN230620A", "HN170419A")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("#Create figure directories")
if(!dir.exists(str_glue('{fig_dir}/{sample_name}'))){
  dir.create(str_glue('{fig_dir}/{sample_name}'))
}

print("#Load data")
print("load data")
HN<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_integrated_object.rds'))

print("#Change the annotation")
print("#Cancer")
HN$numbat_annotation<-ifelse(((as.character(HN$simple_annotation)=="Cancer") | (as.character(HN$simple_annotation)=="Cancer1")| (as.character(HN$simple_annotation)=="Cancer2")| (as.character(HN$simple_annotation)=="Cancer3")| (as.character(HN$simple_annotation)=="Cancer4")| (as.character(HN$simple_annotation)=="Cancer5")| (as.character(HN$simple_annotation)=="Cancer6")), "Malignant", as.character(HN$simple_annotation))

print("#CD4")
HN$numbat_annotation<-ifelse(((as.character(HN$simple_annotation)=="CD4") | (as.character(HN$simple_annotation)=="CD4_1")| (as.character(HN$simple_annotation)=="CD4_2")| (as.character(HN$simple_annotation)=="CD4_3")| (as.character(HN$simple_annotation)=="CD4_4")), "CD4", as.character(HN$numbat_annotation))

print("#CD8")
HN$numbat_annotation<-ifelse(((as.character(HN$simple_annotation)=="CD8") | (as.character(HN$simple_annotation)=="CD8_1")| (as.character(HN$simple_annotation)=="CD8_2")| (as.character(HN$simple_annotation)=="CD8_3")| (as.character(HN$simple_annotation)=="CD8_4")), "CD8", as.character(HN$numbat_annotation))

print("#BCell")
HN$numbat_annotation<-ifelse(((as.character(HN$simple_annotation)=="BCell") | (as.character(HN$simple_annotation)=="BCells")| (as.character(HN$simple_annotation)=="BCell1")| (as.character(HN$simple_annotation)=="BCell2")| (as.character(HN$simple_annotation)=="BCell3")| (as.character(HN$simple_annotation)=="BCell4")), "BCell", as.character(HN$numbat_annotation))

print("#UMAP plot")
Idents(HN)<-"numbat_annotation"
p1<- DimPlot(HN, reduction = "umap", group.by="Location")
p2<-DimPlot(HN, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf(str_glue("{fig_dir}/{sample_name}/numbat_UMAP.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

print("#Save")
saveRDS(HN, file = str_glue("{data_dir}/{sample_name}/{sample_name}_numbat_object.rds"))

