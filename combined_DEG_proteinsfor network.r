#Combine all DEG and do protein network analysis

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


#Everything is named weirdly so just upload separately.

a<-read.csv(file = str_glue("{data_dir}/HN200519A_DEG_analysis_LN_v_primary4.csv"), header = TRUE)
b<-read.csv(file = str_glue("{data_dir}/HN200519A_DEG_analysis_LN_v_primary6.csv"), header = TRUE)
c<-read.csv(file = str_glue("{data_dir}/HN120819A_DEG_analysis_LN2_vP2.csv"), header = TRUE)
d<-read.csv(file = str_glue("{data_dir}/HN120819A_DEG_analysis_LN23_v_P23.csv"), header = TRUE)
e<-read.csv(file = str_glue("{data_dir}/HN021219A_DEG_analysis_LN3_v_primary.csv"), header = TRUE)
f<-read.csv(file = str_glue("{data_dir}/HN021219A_DEG_analysis_LN4_v_primary.csv"), header = TRUE)
g<-read.csv(file = str_glue("{data_dir}/HN230620A_DEG_analysis_LN_v_primary.csv"), header = TRUE)
h<-read.csv(file = str_glue("{data_dir}/HN120819A_DEG_analysis_LN3_v_P3.csv"), header = TRUE)


#Rbind them all

df<-do.call(rbind, list(a,b,c,d,e,f,g,h))

#Filter based on adjusted p.value
df1<-filter(df, p_val_adj<0.05)

#Create vector of list of genes
x<-df1$X

#List unique values
x<-unique(x)

#

#Saveout
write.csv(x, file = str_glue("{data_dir}/combined_DEG_list_for_network.csv"))