#Need to do the clone plots nicely for the paper.

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
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/clone_plots"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

print("set sample names")
samples <- c("HN021219A", "HN120819A", "HN200519A", "HN230620A", "HN070219A", "HN170419A", "HN120520A")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("load in the data")
x<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_clones_final_immune_object.rds'))

print("change clone annotation")
x$cloneannotation<-ifelse(as.character(x$annotation == "1"), "Clone 1", as.character(x$annotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "2"), "Clone 2", as.character(x$cloneannotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "3"), "Clone 3", as.character(x$cloneannotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "4"), "Clone 4", as.character(x$cloneannotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "5"), "Clone 5", as.character(x$cloneannotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "6"), "Clone 6", as.character(x$cloneannotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "7"), "Clone 7", as.character(x$cloneannotation))
x$cloneannotation<-ifelse(as.character(x$annotation == "8"), "Clone 8", as.character(x$cloneannotation))

print("try and plot UMAPs now")
DefaultAssay(x)<-"SCT"
Idents(x)<-"cloneannotation"
x$Location<-factor(x$Location, levels =c("Primary", "Lymph_node"))

print("main combined umap")
pdf(str_glue("{fig_dir}/{sample_name}_combined_UMAP_clones.pdf"), width=11.6, height=8.2)
DimPlot(x, reduction = "umap", split.by = "Location", label = FALSE, pt.size = 0.5) & NoAxes() 
dev.off()

print("combined umap with no titles or legends")
pdf(str_glue("{fig_dir}/{sample_name}_combined_UMAP_clones_nolegend_notitle.pdf"), width=11.6, height=8.2)
p<-DimPlot(x, reduction = "umap", split.by = "Location", label = FALSE, pt.size = 0.5) & NoAxes() & NoLegend() 
p+ theme(strip.text.x =element_text(size = 0))
dev.off()

print("Just do cancer cells, load object ")
y<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds'))


print("change nomenclature")
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone1"), "Clone 1", as.character(y$annotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone2"), "Clone 2", as.character(y$cloneannotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone3"), "Clone 3", as.character(y$cloneannotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone4"), "Clone 4", as.character(y$cloneannotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone5"), "Clone 5", as.character(y$cloneannotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone6"), "Clone 6", as.character(y$cloneannotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone7"), "Clone 7", as.character(y$cloneannotation))
y$cloneannotation<-ifelse(as.character(y$annotation == "CancerClone8"), "Clone 8", as.character(y$cloneannotation))

print("plot cancer clones")
DefaultAssay(y)<-"SCT"
Idents(y)<-"cloneannotation"
y$Location<-factor(y$Location, levels =c("Primary", "Lymph_node"))
pdf(str_glue("{fig_dir}/{sample_name}_cancer_UMAP_clones.pdf"), width=11.6, height=8.2)
DimPlot(y, reduction = "umap", split.by = "Location", label = FALSE, pt.size = 0.5) & NoAxes() 
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_cancer_UMAP_clones_nolegend_notitle.pdf"), width=11.6, height=8.2)
p<-DimPlot(y, reduction = "umap", split.by = "Location", label = FALSE, pt.size = 0.5) & NoAxes() & NoLegend() 
p+ theme(strip.text.x =element_text(size = 0))
dev.off()