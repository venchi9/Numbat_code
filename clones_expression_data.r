#Try and look at cancer clones and see if there are any expression groups which are interesting.

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
cd /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data
R

library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(Nebulosa)
library(glue)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"
output<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"

#Data is here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data

sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"

#Load the data:
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds"))

#Redo the all the SCT, find neibours and clustering again.
x<-SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
x<- RunPCA(x, verbose = FALSE)
#Elbow plot to get dims
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_elbow_canceronly.pdf"), width=11.6, height=8.2)
ElbowPlot(x)
dev.off()
#Cluster
x<-FindNeighbors(x, dims = 1:10)
x<-FindClusters(x, resolution =0.2)
x<-RunUMAP(x, dims = 1:10)



#Do a UMAP
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_UMAP_basic_new.pdf"), width=11.6, height=8.2)
p<-DimPlot(x, reduction = "umap", split.by = "Location")
p
dev.off()

#Need to redo all the SCT, find neibours and clustering again.  Set the resolution at 0.2
Idents(x)<-"clone_opt"
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_UMAP_basic_byclones.pdf"), width=11.6, height=8.2)
p<-DimPlot(x, reduction = "umap", split.by = "Location")
p
dev.off()

#Save object
saveRDS(x, file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#Make the plots more fancy
#Readfile
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
#Need to redo all the SCT, find neibours and clustering again.  Set the resolution at 0.2
Idents(x)<-"clone_opt"
#Reorder based on location
x$Location<-factor(x = x$Location, levels = c("Primary", "Lymph_node"))
pdf(str_glue("{fig_dir}/{sample_name}/{sample_name}_UMAP_basic_byclones_changedcolours.pdf"), width=5.8, height=4.1)
p<-DimPlot(x, reduction = "umap", split.by = "Location", pt.size = 0.8,cols = c("1" = "#2A348D", "2" = "#E41A1C", "3" = "#3F918B", "4" = "#896191", "5" = "#FF980A", "6" ="#70AD47", "7" ="#FFD422"), order = c("7", "6", "5", "4", "3", "2", "1"))
p
dev.off()


#Colours
#1 = #2A348D
#2 = #E41A1C
#3 = #3F918B
#4 = #896191
#5 = #FF980A
#6 = #70AD47
#7 = #FFD422

