#Define the expanding clones in HN200519A, HN021219A, HN120819A, HN230620A
fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate r_4
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

#HN200519A
sample_name<-"HN200519A"

#Load object
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A/HN200519A_clones_final_immune_object.rds")

#Subset the cancer cells only
#First make a meta-data column
hn$cancer<-ifelse((hn$annotation == "1" | hn$annotation == "2" | hn$annotation == "3" | hn$annotation == "4" | hn$annotation == "5" | hn$annotation == "6" | hn$annotation == "7"), "cancer", "non-cancer")

#Now subset
cancer<-subset(hn, subset = cancer == "cancer")

#Write something before the clones
cancer$annotation<-sub("^", "CancerClone", cancer$annotation)

#Add in expanding clones

cancer$expandingclones<-ifelse((cancer$annotation =="CancerClone3" | cancer$annotation =="CancerClone5"), "Expanding", "Non-Expanding")

#Save the cancer object
saveRDS(cancer, file = (str_glue("{data_dir}/{sample_name}_cancerclones_only.rds")))

#Now to the gene lists
#All Ex v non Ex
Idents(cancer)<-"expandingclones"
ex.markers<-FindMarkers(cancer, ident.1 = "Expanding", ident.2 = "Non-Expanding", min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}_ex_v_non_ex.csv")))

#Only the sentinel "nodes"
Idents(cancer)<-"clone_opt"
ex.markers<-FindMarkers(cancer, ident.1 = 3, ident.2 = 6, min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}_ex_sentinel.csv")))

#HN021219A
sample_name<-"HN021219A"

#Load object
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN021219A/HN021219A_clones_final_immune_object.rds")

#Subset the cancer cells only
#First make a meta-data column
hn$cancer<-ifelse((hn$annotation == "1" | hn$annotation == "2" | hn$annotation == "3" | hn$annotation == "4" | hn$annotation == "5" | hn$annotation == "6" | hn$annotation == "7"), "cancer", "non-cancer")

#Now subset
cancer<-subset(hn, subset = cancer == "cancer")

#Write something before the clones
cancer$annotation<-sub("^", "CancerClone", cancer$annotation)

#Add in expanding clones

cancer$expandingclones<-ifelse((cancer$annotation =="CancerClone3" | cancer$annotation =="CancerClone4"), "Expanding", "Non-Expanding")

#Save the cancer object
saveRDS(cancer, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds")))

#Now to the gene lists
#All Ex v non Ex
Idents(cancer)<-"expandingclones"
ex.markers<-FindMarkers(cancer, ident.1 = "Expanding", ident.2 = "Non-Expanding", min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_ex_v_non_ex.csv")))

#Only the sentinel "nodes"
Idents(cancer)<-"clone_opt"
ex.markers<-FindMarkers(cancer, ident.1 = c(3, 4), ident.2 = c(1, 2, 5, 6, 7), min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_ex_sentinel.csv")))

#HN120819A
sample_name<-"HN120819A"

#Load object
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/HN120819A_clones_final_immune_object.rds")

#Subset the cancer cells only
#First make a meta-data column
hn$cancer<-ifelse((hn$annotation == "1" | hn$annotation == "2" | hn$annotation == "3" | hn$annotation == "4" | hn$annotation == "5" | hn$annotation == "6" | hn$annotation == "7"), "cancer", "non-cancer")

#Now subset
cancer<-subset(hn, subset = cancer == "cancer")

#Write something before the clones
cancer$annotation<-sub("^", "CancerClone", cancer$annotation)

#Add in expanding clones

cancer$expandingclones<-ifelse((cancer$annotation =="CancerClone3"), "Expanding", "Non-Expanding")

#Save the cancer object
saveRDS(cancer, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds")))

#Now to the gene lists
#All Ex v non Ex
Idents(cancer)<-"expandingclones"
ex.markers<-FindMarkers(cancer, ident.1 = "Expanding", ident.2 = "Non-Expanding", min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_ex_v_non_ex.csv")))

#Only the sentinel "nodes"
Idents(cancer)<-"clone_opt"
ex.markers<-FindMarkers(cancer, ident.1 = c(3), ident.2 = c(1, 2), min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_ex_sentinel.csv")))

#HN230620A
sample_name<-"HN230620A"

#Load object
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN230620A/HN230620A_clones_final_immune_object.rds")

#Subset the cancer cells only
#First make a meta-data column
hn$cancer<-ifelse((hn$annotation == "1" | hn$annotation == "2" | hn$annotation == "3" | hn$annotation == "4" | hn$annotation == "5" | hn$annotation == "6" | hn$annotation == "7"), "cancer", "non-cancer")

#Now subset
cancer<-subset(hn, subset = cancer == "cancer")

#Write something before the clones
cancer$annotation<-sub("^", "CancerClone", cancer$annotation)

#Add in expanding clones

cancer$expandingclones<-ifelse((cancer$annotation =="CancerClone3"), "Expanding", "Non-Expanding")

#Save the cancer object
saveRDS(cancer, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_cancerclones_only.rds")))

#Now to the gene lists
#All Ex v non Ex
Idents(cancer)<-"expandingclones"
ex.markers<-FindMarkers(cancer, ident.1 = "Expanding", ident.2 = "Non-Expanding", min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_ex_v_non_ex.csv")))

#Only the sentinel "nodes"
Idents(cancer)<-"clone_opt"
ex.markers<-FindMarkers(cancer, ident.1 = c(4, 5), ident.2 = c(1, 2, 3, 6), min.pct = 0.25, vars.to.regress = "orig.ident")
#write
write.csv(ex.markers, file = (str_glue("{data_dir}/{sample_name}/{sample_name}_ex_sentinel.csv")))




