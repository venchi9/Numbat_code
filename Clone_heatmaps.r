#Want to look at heatmaps for each clone to work out highly differential genes as a start


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

print("load directories")
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"
output<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"

print("load in samples")
samples <- c("HN230620A",  "HN170419A", "HN120819A", "HN120520A", "HN070219A", "HN021219A") #I need to do "HN200519A" later as the clones aren't finalised yet.
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("#Create figure directories")
if(!dir.exists(str_glue('{fig_dir}/{sample_name}'))){
  dir.create(str_glue('{fig_dir}/{sample_name}'))
}

print("read in object")
hn<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_numbat_object_noduplicates_clones.rds'))

print("create object of only cancer cells")
cancer<-subset(hn, subset = clone_opt != "NA")

print("find variable markers")
Idents(cancer)<-"clone_opt"
DefaultAssay(cancer)<-"RNA"
print("scale data")
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
print("find variable markers")
pbmc.markers <- FindAllMarkers(cancer, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

write.csv(pbmc.markers, str_glue("{fig_dir}/{sample_name}/top10markers.csv"))

print("draw heatmap")
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf(str_glue("{fig_dir}/{sample_name}/clone_heatmap.pdf"), width=11.6, height=8.2)
DoHeatmap(cancer, features = top10$gene) + NoLegend() +theme(text=element_text(size=6)) 
dev.off()


