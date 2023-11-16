#DEG analysis of the transcriptional clonal groups
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

#Start with HN200519A
sample_name<-"HN200519A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
#Two comparisons - lymph node expanding versus primary clone 4 and primary clone 6

#Comparison 1 - LN expanding v primary clone 4
Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_expanding", ident.2 = "Primary_4", vars.to.regress = orig.ident)

#Save out
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN_v_primary4.csv"))

#Comparison 2 - LN expanding v primary clone 6

Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_expanding", ident.2 = "Primary_6", vars.to.regress = orig.ident)

#Save out
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN_v_primary6.csv"))

#HN120819A
#For this group we want primary 23 v LN 23 and probably P2 v LN2 and P3 v LN3
#comparison 1:  P2 v LN2
sample_name<-"HN120819A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_2", ident.2 = "primary_2", vars.to.regress = orig.ident)

#Save out
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN2_vP2.csv"))
#Comparison 2:  LN3 v P3
Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_3", ident.2 = "primary_3", vars.to.regress = orig.ident)
#Save out
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN3_v_P3.csv"))

#Comparison 3 - combined the LN and primary groups
x$comparison<-ifelse(((x$DEGcluster =="Lymph_2") | (x$DEGcluster =="Lymph_3")), "LN", as.character(x$DEGcluster))
x$comparison<-ifelse(((x$DEGcluster =="primary_2") | (x$DEGcluster =="primary_3")), "primary", as.character(x$comparison))

Idents(x)<-"comparison"
exmarkers<-FindMarkers(x, ident.1 = "LN", ident.2 = "primary", vars.to.regress = orig.ident)

#saveout
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN23_v_P23.csv"))

#HN021219A
#Comparisons 1: are clone 3 LN v primary 5,6,7 and clone 4 LN v primary 5,6,7

sample_name<-"HN021219A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))

#Comparison 1:  primary v Lymph 3
Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_3", ident.2 = "primary", vars.to.regress = orig.ident)

#saveout
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN3_v_primary.csv"))

#Comparison 2:  lymph 4 v primary
Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_4", ident.2 = "primary", vars.to.regress = orig.ident)
#saveout
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN4_v_primary.csv"))

#HN230620A
sample_name<-"HN230620A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
#Only one comparison - clone 45 (LN) v clone 2,3,6 (primary)
Idents(x)<-"DEGcluster"
exmarkers<-FindMarkers(x, ident.1 = "Lymph_node", ident.2 = "primary", vars.to.regress = orig.ident)
#saveout
write.csv(exmarkers, file = str_glue("{data_dir}/{sample_name}_DEG_analysis_LN_v_primary.csv"))



#Now do the combined object which I painfully curated

#object is here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/merged_object_paper_SCT.rds

Idents(merged)<-"final"
exmarkers<-FindMarkers(merged, ident.1 = "LN", ident.2 = "Primary", vars.to.regress = orig.ident)
ordered<-exmarkers[order(exmarkers$avg_log2FC),]

write.csv(exmarkers,file = str_glue("{data_dir}/combined_cancer_DEG.csv"))