#Take the final cancer only objects and create DEG lists for the expanding versus non-expanding clones.
#Here is what I've agreed on with Joseph as the comparisons:

#HN200519A - 3 + 5 v everything else
#HN120819A - 3 v everything else
#HN021219A - 4 + 3 v everything else
#HN230620A - 4+5 v everything else.

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2

library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(Nebulosa)
library(glue)
library(MAST)
library(DESeq2)

#Data is here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data
#Can put figures here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures
#Gene lists here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
gene_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists"

#Define sample name
sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"
sample_name<-"combined"

#Load in data
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds"))

#Do differential gene expression analysis
Idents(x)<-"expandingclones"
exmarkers<-FindMarkers(x, ident.1 = "Expanding", ident.2 = "Non-Expanding", vars.to.regress = orig.ident)

#Write the data
write.csv(exmarkers, file = str_glue("{gene_dir}/{sample_name}_DEG_list.csv"))

#Done

#Try MAST
gene_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists/MAST"

sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"
sample_name<-"combined"

#Load in data
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds"))

#Do differential gene expression analysis
Idents(x)<-"expandingclones"
DefaultAssay(x)<-"RNA"
x<-PrepSCTFindMarkers(object = x)
exmarkers<-FindMarkers(x, ident.1 = "Expanding", ident.2 = "Non-Expanding", test.use = "t", vars.to.regress = orig.ident)
#There is a bug with using SCT

#Write the data
write.csv(exmarkers, file = str_glue("{gene_dir}/{sample_name}_DEG_list_MAST.csv"))




