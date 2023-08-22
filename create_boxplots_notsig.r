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
samples <- c("EIF4A1", "BAK1", "TAP2")
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
pdf(str_glue("{fig_dir}/boxplot_notsig_{sample_name}_witht_test.pdf"), width=4, height=4)
level_order<-c('Primary', 'LN')
my_comparisons<-list(c("Primary", "LN"))
p <- ggplot(combined, aes( x = factor(Group, level = level_order), y = expression)) +
  geom_boxplot(aes(fill = Group)) 
 p + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + theme_bw()+ scale_fill_manual(values = c("palegreen1", "green4"))
 dev.off()

print("plot without legends")
 pdf(str_glue("{fig_dir}/boxplot_notsig_{sample_name}_nolegends.pdf"), width=3, height=3)
level_order<-c('Primary', 'LN')
my_comparisons<-list(c("Primary", "LN"))
p <- ggplot(combined, aes( x = factor(Group, level = level_order), y = expression)) +
  geom_boxplot(aes(fill = Group)) + ylim(-5, 6)
 p  + theme_bw() + NoLegend() + theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+ scale_fill_manual(values = c("palegreen1", "green4"))
 dev.off()