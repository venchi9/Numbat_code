#Try and create a VENN diagram for the CNV areas between individuals
#Data is here /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data
#Individuals I'm interested in:  HN021219A, HN120819A, HN200519A, HN230620A
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

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures/clone_plots"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

#First load up the data:
#HN021219A
HN021219A<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN021219A/numbat/segs_consensus_20.tsv")
HN021219A_list<-HN021219A$seg
HN021219A_list<-HN021219A_list[!is.na(HN021219A_list)]



#HN120819A
HN120819A<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN120819A/numbat/segs_consensus_20.tsv")
HN120819A_list<-HN120819A$seg
HN120819A_list<-HN120819A_list[!is.na(HN120819A_list)]


#HN200519A
HN200519A<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN200519A/numbat_400/segs_consensus_20.tsv")
HN200519A_list<-HN200519A$seg
HN200519A_list<-HN200519A_list[!is.na(HN200519A_list)]


#HN230620A
HN230620A<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN230620A/numbat/segs_consensus_20.tsv")
HN230620A_list<-HN230620A$seg
HN230620A_list<-HN230620A_list[!is.na(HN230620A_list)]

#Need to make all vectors the same length
#Longest is HN021219A (n=25)
n<-max(length(HN021219A_list))
length(HN200519A_list)<-n
length(HN230620A_list)<-n
length(HN120819A_list)<-n
#Change to data frame
HN021219A_list<-as.data.frame(HN021219A_list)
HN200519A_list<-as.data.frame(HN200519A_list)
HN230620A_list<-as.data.frame(HN230620A_list)
HN120819A_list<-as.data.frame(HN120819A_list)
#Add metadata column to include individual name
HN021219A_list$individual<-"HN021219A"
HN200519A_list$individual<-"HN200519A"
HN230620A_list$individual<-"HN230620A"
HN120819A_list$individual<-"HN120819A"
#Make the column names identical
colnames(HN021219A_list)<-c("segments", "individual")
colnames(HN200519A_list)<-c("segments", "individual")
colnames(HN230620A_list)<-c("segments", "individual")
colnames(HN120819A_list)<-c("segments", "individual")
#Create a single table
#Bind the dataframes together
CNV<-rbind(HN021219A_list, HN120819A_list, HN200519A_list, HN230620A_list)

#Let's try and pivot this table
test<-CNV%>%
    select(CNV$individual, CNV$segments) %>%
    summarise(count=n()) %>%
    pivot_wider(names_from = CNV$segments,
                values_from = count)%>%
    replace(is.na(.), 0)%>%
    as.data.frame()


> HN021219A_list<-HN021219A$seg
> HN120819A_list<-HN120819A$seg
> HN200519A_list<-HN200519A$seg
> HN230620A_list<-HN230620A$seg
> listInput<-list(patient_1 =  HN200519A_list, patient_2 = HN120819A_list, patient_3 = HN021219A_list,   patient_4 =  HN230620A_list)
pdf(str_glue("{fig_dir}/upsetR_plot.pdf"), width=11.6, height=8.2)
upset(fromList(listInput), order.by = "freq")
dev.off()

#Make into matrix
library(qdapTools)
listInput<-list(patient_1 =  HN200519A_list, patient_2 = HN120819A_list, patient_3 = HN021219A_list,   patient_4 =  HN230620A_list)
mtabulate(listInput)
listInput<-mtabulate(listInput)

#Plot
pdf(str_glue("{fig_dir}/upsetR_plot.pdf"), width=11.6, height=8.2)
upset(listInput, sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")
dev.off()
#That works.

#Try now with some better labelling

#Try a VEN diagram
library("ggVennDiagram")
listInput<-list(patient_1 =  HN200519A_list, patient_2 = HN120819A_list, patient_3 = HN021219A_list,   patient_4 =  HN230620A_list)
pdf(str_glue("{fig_dir}/VENN_diagram.pdf"), width=11.6, height=8.2)
ggVennDiagram(listInput)
dev.off()
#Try with some extra fancy things
pdf(str_glue("{fig_dir}/VENN_diagram_fancy.pdf"), width=11.6, height=8.2)
ggVennDiagram(
  listInput, label_alpha = 0,
  category.names = c("patient_1","patient_2","patient_3", "patient_4")
  ) +
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")
  dev.off()

  #Try ggvenn because I can't load up the one above
  pdf(str_glue("{fig_dir}/ggVENN_diagram.pdf"), width=11.6, height=8.2)
  listInput<-list(patient_1 =  HN200519A_list, patient_2 = HN120819A_list, patient_3 = HN021219A_list,   patient_4 =  HN230620A_list)
  ggvenn(listInput)
  dev.off()

  #Find the intersection of the gene lists
a<-intersect(HN200519A_list, HN120819A_list)
#"6a"  "8a"  "16a"
b<-intersect(HN200519A_list, HN021219A_list)
#"6a"  "8a"  "16c"
c<-intersect(HN200519A_list, HN230620A_list)
#"20a"
d<-intersect(HN120819A_list, HN200519A_list)
#"6a"  "8a"  "16a"
e<-intersect(HN120819A_list, HN021219A_list)
#"2a"  "6a"  "8a"  "11a" "11d" "12b"
f<-intersect(HN120819A_list, HN230620A_list)
#"8b"
g<-intersect(HN021219A_list, HN200519A_list)
#"6a"  "8a"  "16c"
h<-intersect(HN021219A_list, HN120819A_list)
#"2a"  "6a"  "8a"  "11a" "11d" "12b"
i<-intersect(HN021219A_list, HN230620A_list)
#"3b" "5a" "7b"
j<-intersect(HN230620A_list, HN200519A_list)
#"20a"
k<-intersect(HN230620A_list, HN120819A_list)
#"8b"
l<-intersect(HN230620A_list, HN021219A_list)
#"3b" "5a" "7b"

#Next, generate a list of all the patients, with the CNVs and what type they are.  





