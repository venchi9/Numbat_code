#Creating boxplots for all the genes I'm interested in

suppressPackageStartupMessages({
library(tidyverse)
library(cowplot)
library(tidyverse)
library(glue)
library(Seurat)
library(scCustomize)
library(stats)
library(rstatix)
library(ggpubr)
})

#Load in gene names
print("load in samples")
samples <- c("EIF4G1",   "EIF4A1" ,  "EIF4E"  ,  "EIF4EBP1" ,"PDCD4"  ,  "EIF2S1"  ,
"EIF2AK1",  "EIF2AK2" , "EIF2AK3" , "EIF2AK4" , "PPP1R15A", "PPP1CA" , 
"VEGFA" ,   "HIF1A" ,   "FGF2"  ,   "BCL2" ,    "TRPV6"  ,  "SLC38A2" ,
"XIAP" ,    "SREBF1" ,  "BIRCR"  ,  "GULP"  ,   "PCBP1" ,   "PCBP2"   ,
 "EIF4G2",   "EIF4G3"  , "EIF4A2" ,  "EIF4A3"  , "MTOR"  ,   "ATP7A"   ,
"EEF2"  ,   "RPS6"  ,   "EEF2K"  ,  "MYC"  ,    "TP53"  ,  "RB1"    , 
"RPL36A" ,  "RPL23A"  , "RL35A" ,   "RPL10"  ,  "RPL5" ,    "RPL39"  , 
"RPL22"  ,  "RPS15" ,   "RPS17"  ,  "RPS24"  ,  "RPS19" ,   "HLA-C"  , 
"HLA-B"  ,  "HLA-A" ,   "B2M"  ,    "TAP1" ,    "TAP2" ,    "TAPBP"  , 
"HLA-DRB1", "HLA-DRA" , "DLA-DRB3", "HLA-DRB4" ,"HLA-DQB1" ,"HLA-DQA1", "RHOA", "SLIT1", "ROBO2", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10")
#n = 60
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]
#Set figure directory
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures/geneboxplots"

print("Read in RDS file")
allmerged<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined_cancer_allmerged_allfeatures.rds")

print("extract the gene matrix")

groupA<-subset(allmerged, subset = final =="Primary")
groupB<-subset(allmerged, subset = final =="LN")

expression_A<-data.frame(expression = groupA@assays$SCT@scale.data[sample_name,])
expression_B<-data.frame(expression = groupB@assays$SCT@scale.data[sample_name,])

print("combine the data")
#Combine the data
expression_A$Group<-"Primary"
expression_B$Group<-"LN"
combined<-rbind(expression_A, expression_B)

print("plot with legends")
pdf(str_glue("{fig_dir}/boxplot_{sample_name}_witht_test.pdf"), width=4, height=4)
level_order<-c('Primary', 'LN')
my_comparisons<-list(c("Primary", "LN"))
p <- ggplot(combined, aes( x = factor(Group, level = level_order), y = expression)) +
  geom_boxplot(aes(fill = Group)) 
 p + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + theme_bw()
 dev.off()

print("plot without legends")
 pdf(str_glue("{fig_dir}/boxplot_{sample_name}_nolegends.pdf"), width=3, height=3)
level_order<-c('Primary', 'LN')
my_comparisons<-list(c("Primary", "LN"))
p <- ggplot(combined, aes( x = factor(Group, level = level_order), y = expression)) +
  geom_boxplot(aes(fill = Group)) + ylim(-5, 6)
 p  + theme_bw() + NoLegend() + theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
 dev.off()
