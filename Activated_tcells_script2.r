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
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/code"

samples <- c("LN14_filter_SCT_PCA_simpleannotation", "LN23_rep1_filter_SCT_PCA_simpleannotation", "LN23_rep2_filter_SCT_PCA_simpleannotation", "pri34_filter_SCT_PCA_simpleannotation", "pri234_filter_SCT_PCA_simpleannotation", "RZ726_LN2_filter_SCT_PCA_simpleannotation", "RZLN_filter_SCT_PCA_simpleannotation","RZpri_filter_SCT_PCA_simpleannotation")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("read in object")
immune<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_macrophage_immune_azimuth_annotation_tcell.rds'))
#Do a UMAP
p1<- DimPlot(immune, reduction = "umap", group.by = "MajoritySinglet_Individual_Assignment", label = TRUE, repel = TRUE) 
p2<-DimPlot(immune, reduction = "umap", group.by = "final_immune", label = TRUE, repel = TRUE) + NoLegend()
pdf(str_glue("{fig_dir}/{sample_name}/UMAP_final_immune_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Subset out the individuals and then create .csv files for them to put back in master objects
table(immune$MajoritySinglet_Individual_Assignment)
p1<-subset(immune, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
p2<-subset(immune, subset = MajoritySinglet_Individual_Assignment == "HN200519A")
p3<-subset(immune, subset = MajoritySinglet_Individual_Assignment == "HN230620A")
#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)

annotation<-p3$final_immune
df3<-data.frame(annotation)

#Save out
Ind_name<-"HN170419A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN200519A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN230620A"
write.csv(df3, file = str_glue("{data_dir}/{sample_name}_{Ind_name}__immunecells_annotation.csv"))

#Done.  Now just need to insert into each big object.  That will be very fun.

#Okay, let's start with HN021219A
#The main object is here:
#/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN021219A/HN021219A_numbat_object_noduplicates_clones.rds
#Immune cell labeling from:
#RZ726_LN2:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_LN2_HN021219A_immunecells_annotation.csv
#RZ726_pri:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_pri_HN021219A_immunecells_annotation.csv

sample_name<-"HN021219A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN021219A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN021219A"

#Read in the object:
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN021219A/HN021219A_numbat_object_noduplicates_clones.rds")
immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_pri_HN021219A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_LN2_HN021219A_immunecells_annotation.csv")
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X

#bind the .csv files togehter:

comb<-rbind(immune1, immune2)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))
#worked!!!!

#HN070219A
#Main object:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN070219A/HN070219A_numbat_object_noduplicates_clones.rds

sample_name<-"HN070219A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN070219A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN070219A"

#Immune cell labeling from:
#LN23_rep1:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep1_filter_SCT_PCA_simpleannotation_HN070219A_immunecells_annotation.csv
#LN23_rep2:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep2_filter_SCT_PCA_simpleannotation_HN070219A_immunecells_annotation.csv

hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN070219A/HN070219A_numbat_object_noduplicates_clones.rds")
immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep1_HN070219A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep2_HN070219A_immunecells_annotation.csv")

#Assign rownames
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X

#bind the .csv files togehter:

comb<-rbind(immune1, immune2)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Get rid of excess columns
hn$X<-NULL
hn$cell<-NULL

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))


#HN120520A
#Main object:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120520A/HN120520A_numbat_object_noduplicates_clones.rds

sample_name<-"HN120520A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120520A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN120520A"


#Immune cell labeling from:
#RZ726_pri: /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZpri_filter_SCT_PCA_simpleannotation_HN120520A_immunecells_annotation.csv
#RZ726_LN:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZLN_filter_SCT_PCA_simpleannotation_HN120520A_immunecells_annotation.csv


hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120520A/HN120520A_numbat_object_noduplicates_clones.rds")
immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_pri_HN120520A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_LN_HN120520A_immunecells_annotation.csv")

#Assign rownames
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X

#bind the .csv files togehter:

comb<-rbind(immune1, immune2)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Get rid of excess columns
hn$X<-NULL
hn$cell<-NULL

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))

#HN120819A
#Main object:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/HN120819A_numbat_object_noduplicates_clones.rds

sample_name<-"HN120819A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN120819A"

#Immune cell labeleing from:
#pri_234:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri234_filter_SCT_PCA_simpleannotation_HN120819A_immunecells_annotation.csv
#pri_34:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri34_filter_SCT_PCA_simpleannotation_HN120819A_immunecells_annotation.csv
#LN23_rep1:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep1_filter_SCT_PCA_simpleannotation_HN120819A_immunecells_annotation.csv
#LN23_rep2:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep2_filter_SCT_PCA_simpleannotation_HN120819A_immunecells_annotation.csv

hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/HN120819A_numbat_object_noduplicates_clones.rds")
immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_234_HN120819A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_34_HN120819A_immunecells_annotation.csv")
immune3<-read.csv(file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep1_HN120819A_immunecells_annotation.csv")
immune4<-read.csv(file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep2_HN120819A_immunecells_annotation.csv")

#Assign rownames
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X
rownames(immune3)<-immune3$X
rownames(immune4)<-immune4$X

#bind the .csv files togehter:

comb<-rbind(immune1, immune2, immune3, immune4)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Get rid of excess columns
hn$X<-NULL
hn$cell<-NULL

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))

#HN170419A
#main object:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN170419A/HN170419A_numbat_object_noduplicates_clones.rds

sample_name<-"HN170419A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN170419A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN170419A"

#Immune cell labeling from:
#pri_234:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri234_filter_SCT_PCA_simpleannotation_HN170419A_immunecells_annotation.csv
#LN14:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_filter_SCT_PCA_simpleannotation_HN170419A_immunecells_annotation.csv

hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN170419A/HN170419A_numbat_object_noduplicates_clones.rds")
immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_234_HN170419A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_HN170419A_immunecells_annotation.csv")

#Assign rownames
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X

comb<-rbind(immune1, immune2)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Get rid of excess columns
hn$X<-NULL
hn$cell<-NULL

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))

#HN230620A

sample_name<-"HN230620A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN230620A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN230620A"

#Main object:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN230620A/HN230620A_numbat_object_noduplicates_clones.rds
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN230620A/HN230620A_numbat_object_noduplicates_clones.rds")

#Immune labeling from:
#RZ726_pri:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZpri_filter_SCT_PCA_simpleannotation_HN230620A__immunecells_annotation.csv
#RZ726_LN:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZLN_filter_SCT_PCA_simpleannotation_HN230620A_immunecells_annotation.csv

immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_pri_HN230620A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ_LN_HN230620A_immunecells_annotation.csv")

#Assign rownames
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X

comb<-rbind(immune1, immune2)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Get rid of excess columns
hn$X<-NULL
hn$cell<-NULL

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))


#HN200519A - need to wait until clones are done.
sample_name<-"HN200519A"
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/HN200519A"

#Main object:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A/HN200519A_numbat400_object_noduplicates_clones.rds
hn<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A/HN200519A_numbat400_object_noduplicates_clones.rds")

#Immune labeling from:
#pri_234:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_234_HN200519A_immunecells_annotation.csv
#pri_34:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_34_HN200519A_immunecells_annotation.csv
#LN14: /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_HN200519A_immunecells_annotation.csv

immune1<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_234_HN200519A_immunecells_annotation.csv")
immune2<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri_34_HN200519A_immunecells_annotation.csv")
immune3<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_HN200519A_immunecells_annotation.csv")
#Assign rownames
rownames(immune1)<-immune1$X
rownames(immune2)<-immune2$X
rownames(immune3)<-immune3$X

comb<-rbind(immune1, immune2, immune3)
hn<-AddMetaData(object = hn, metadata = comb)

#Put in the clones

hn$annotation<-ifelse((is.na(hn$annotation)), hn$clone_opt, as.character(hn$annotation))
#put in stroma and endothelial cells
hn$annotation<-ifelse((is.na(hn$annotation)), hn$numbat_annotation, as.character(hn$annotation))

#Get rid of excess columns
hn$X<-NULL
hn$cell<-NULL

#Do a UMAP
Idents(hn)<-"annotation"
p1<- DimPlot(hn, reduction = "umap", group.by="annotation", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(hn, reduction = "umap", group.by = "Location") 
pdf(str_glue("{fig_dir}/{sample_name}_final_UMAP_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(hn, file = (str_glue("{data_dir}/{sample_name}_clones_final_immune_object.rds")))



test1<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_filter_SCT_PCA_simpleannotation/LN14_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")
test2<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_filter_SCT_PCA_simpleannotation/LN14_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")



