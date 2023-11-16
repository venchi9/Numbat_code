#Find which of the curated genes are in the enriched pathways
fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(cowplot)
library(tidyverse)
library(glue)
library(Seurat)
library(scCustomize)
library(stats)
library(rstatix)
library(ggpubr)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
list_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/reactomeGeneLists"


#The top enriched pathways
#ViralmRNA
pathway<-read_tsv(file = str_glue("{data_dir}/viralmRNA_genelist.txt"), col_names = FALSE)
pathway<-pathway$X1
#protein translation list
protein_list<-c("TP53", "BAX", "FAS", "BAK1","UBE3A", "EIF2S1","EIF2AK2","PPP1CC","PPP1CB","PPP1CA","PPP1R15A")

#Find intersection
a<-intersect(pathway, protein_list)
#Nil

#Immune evasion list
gene_list<-c("TNF", "LTA", "LTB", "BAK1", "CFB", "STK19", "HSPA1A", "HSP90AA1", "HSP90B1", "HLA-A", "HLA-B", "HLA-C", "IRF3",  "IRF7", "MX2", "IFITM3", "ISG15","OASL", "JAK1", "JAK2", "STAT1", "STAT2", "TYK2", "IFNAR1", "IFNAR2", "IFNGR1","IFNGR2" ,"IFIH1","NFKB1")
b<-intersect(pathway, gene_list)

#Nil

#Peptide chain
pathway<-read_tsv(file = str_glue("{data_dir}/peptide_chain.txt"), col_names = FALSE)
pathway<-pathway$X1

#Nil intersection

#SRP dependent
pathway<-read_tsv(file = str_glue("{data_dir}/SRP_dependent.txt"), col_names = FALSE)
pathway<-pathway$X1

#Nil intersection


#Regulation SLIT ROBO
pathway<-read_tsv(file = str_glue("{data_dir}/Regulation_slit_robo.txt"), col_names = FALSE)
pathway<-pathway$X1

#Nil

#ROBO
pathway<-read_tsv(file = str_glue("{data_dir}/ROBO_signaling.txt"), col_names = FALSE)
pathway<-pathway$X1
#Nil

#Eukaryotic
pathway<-read_tsv(file = str_glue("{data_dir}/Eukaryotic.txt"), col_names = FALSE)
pathway<-pathway$X1
#Nil

#40S subunit
pathway<-read_tsv(file = str_glue("{data_dir}/40S_subunit.txt"), col_names = FALSE)
pathway<-pathway$X1

#No interaction between genes at all.



