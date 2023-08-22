#Take all the common pathways found in Reactome and get genelists for each pathway for each comparison

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(cowplot)
library(tidyverse)
library(glue)


#Data is here:
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
list_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/reactomeGeneLists"

#The commonly enriched pathways are:

#Viral mRNA Translation
#Signaling by ROBO receptors
#Regulation of expression of SLITs and ROBOs
#Formation of a pool of free 40S subunits
#Peptide chain elongation

#Samples
sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"

#Comparisons
comparison_name<-"LN_v_primary4"
comparison_name<-"LN_v_primary6" - no viral mRNA translation
comparison_name<-"LN2_vP2"
comparison_name<-"LN3_v_P3"
comparison_name<-"LN23_v_P23"
comparison_name<-"LN3_v_primary"
comparison_name<-"LN4_v_primary"
comparison_name<-"LN_v_primary"

#Load up the reactome pathways.
geneset<-read.csv(file = str_glue("{data_dir}/{sample_name}_{comparison_name}_reactome_enriched_pathways.csv"))

#Make vectors for the top 10 genelists
x<-geneset[1:10, ]

#Walter's solution
gene.sets <- strsplit(as.character(x$geneID), "/")
names(gene.sets) <- x$Description

#Find the pathway list
names(gene.sets)

#viral translation
a<-unlist(gene.sets[[7]])


#Signaling ROBO
b<-unlist(gene.sets[[10]])


#Regulation of slit/robo
c<-unlist(gene.sets[[8]])


#Formation of 40S subunit
d<-unlist(gene.sets[[4]])


#Peptide elongation
e<-unlist(gene.sets[[1]])


#Make all vectors the same length
n<-max(length(a),length(b), length(c) ,length(d), length(e))

length(a)<-n
length(b)<-n
length(c)<-n
length(d)<-n
length(e)<-n
#Bind into a single table
df<-cbind(a,b,c,d,e)

#Bind into single data.frame


df<-as.data.frame(df)

#Set column names
setnames(df, old = c("a", "d", "e"), 
         new = c("Viral mRNA Translation","Signaling by ROBO receptors","Regulation of expression of SLITs and ROBOs", "Formation of a pool of free 40S subunits", "Peptide chain elongation"))

#Save out
write.csv(df, file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))