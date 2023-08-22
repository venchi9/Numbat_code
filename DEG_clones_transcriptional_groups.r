#Cancer clones have been re-grouped into transcriptional clusters.  Joseph and I have worked out the comparisons that need to be done
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

#load the sample

#Start with HN200519A

sample_name<-"HN200519A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#First create the groups for comparison.
#For this patient we want two comparisons
#The lymph node group is seurat cluster 0+2+4 and clones 3 + 5
#Primary group one is seruat clusters 3 + 4 and clone 4
#Primary group two is seurat cluster 1 and clone 6

#Create LN group
x$DEGcluster<-ifelse(((x$Location =="Lymph_node")&((x$seurat_clusters == "0")| (x$seurat_clusters == "2") | (x$seurat_clusters == "4")) & ((x$clone_opt == "3") | x$clone_opt == "5")), "Lymph_expanding", "other")

#Create primary group 1
x$DEGcluster<-ifelse(((x$Location =="Primary")&((x$seurat_clusters == "3")| (x$seurat_clusters == "4")) & (x$clone_opt == "4")), "Primary_4", as.character(x$DEGcluster))

#Creat primary group2
x$DEGcluster<-ifelse(((x$Location =="Primary")&(x$seurat_clusters == "1") & (x$clone_opt == "6")), "Primary_6", as.character(x$DEGcluster))

#Save
saveRDS(x, file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#Now HN120819A
sample_name<-"HN120819A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#For this group we want primary 23 v LN 23 and probably P2 v LN2 and P3 v LN3
#Create LN group
x$DEGcluster<-ifelse(((x$Location =="Lymph_node")&((x$seurat_clusters == "3")| (x$seurat_clusters == "0") | (x$seurat_clusters == "2")| (x$seurat_clusters == "4")) & (x$clone_opt == "2")), "Lymph_2", "other")
#Create next LN group
x$DEGcluster<-ifelse(((x$Location =="Lymph_node")&((x$seurat_clusters == "3")| (x$seurat_clusters == "0") | (x$seurat_clusters == "2")| (x$seurat_clusters == "4")) & (x$clone_opt == "3")), "Lymph_3", as.character(x$DEGcluster))
#Create primary group 1
x$DEGcluster<-ifelse(((x$Location =="Primary")&((x$seurat_clusters == "1")| (x$seurat_clusters == "2")) & (x$clone_opt == "2")), "primary_2", as.character(x$DEGcluster))
#Create primary group2
x$DEGcluster<-ifelse(((x$Location =="Primary")&((x$seurat_clusters == "1")| (x$seurat_clusters == "2")) & (x$clone_opt == "3")), "primary_3", as.character(x$DEGcluster))

#save out
saveRDS(x, file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))


#HN021219A
#We want LN group with clone 3
#LN group with clone 4
#Primary group with clone 5, 6 7
sample_name<-"HN021219A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#Create LN group
x$DEGcluster<-ifelse(((x$Location =="Lymph_node")&((x$seurat_clusters == "0")| (x$seurat_clusters == "3")) & (x$clone_opt == "3")), "Lymph_3", "other")
#Second LN group
x$DEGcluster<-ifelse(((x$Location =="Lymph_node")&(x$seurat_clusters == "1") & (x$clone_opt == "4")), "Lymph_4", as.character(x$DEGcluster))
#Primary group
x$DEGcluster<-ifelse(((x$Location =="Primary")&(x$seurat_clusters == "2") & ((x$clone_opt == "5") | (x$clone_opt == "6")|(x$clone_opt == "7"))), "primary", as.character(x$DEGcluster))

#save out
saveRDS(x, file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#HN230620A
sample_name<-"HN230620A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
#Less well defined group
#LN group is seurat cluster 0,1,2 and clones 4,5
#primary is suerat 3, clone 2,3,6

#First primary group
x$DEGcluster<-ifelse(((x$Location =="Primary")&(x$seurat_clusters == "3") & ((x$clone_opt == "2") | (x$clone_opt == "3")| (x$clone_opt == "6"))), "primary", "other")
#LN group
x$DEGcluster<-ifelse(((x$Location =="Lymph_node")&((x$seurat_clusters == "0")| (x$seurat_clusters == "1")| (x$seurat_clusters == "2"))& ((x$clone_opt == "4") | (x$clone_opt == "5"))), "Lymph_node", as.character(x$DEGcluster))

#save out
saveRDS(x, file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#All groups have been made.  Now need to do DEG.  A job for post-Bali brain.
