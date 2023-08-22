#Use intersecting gene lists and do dotplots for each patient and clone
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

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
list_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/reactomeGeneLists"


#Data is here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data

sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"


#Start with HN200519A
sample_name<-"HN200519A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
#Comparison 1 - LN expanding v primary clone 4
#Comparison 2 - LN expanding v primary clone 6
Idents(x)<-"DEGcluster"

#Get rid of "other cells"
x<-subset(x, idents = c("Lymph_expanding", "Primary_4", "Primary_6"))

#HN120819A
#comparison 1:  P2 v LN2
#Comparison 2:  LN3 v P3
#Comparison 3 - combined the LN and primary groups
x<-subset(x, idents = c("Lymph_2", "Lymph_3", "primary_2", "primary_3"))

#HN021219A
#Comparison 1:  primary v Lymph 3
#Comparison 2:  lymph 4 v primary
x<-subset(x, idents = c("Lymph_3", "Lymph_4", "primary"))

#HN230620A
x<-subset(x, idents = c("Lymph_node", "primary"))


#Start with viral list
list_name<-"viralmRNA"
gene_list<-read.csv(file = str_glue("{list_dir}/viralmRNA_intersection_list.csv"))
gene_list<-gene_list$x
gene_list<-na.omit(gene_list)

#Create dotplot
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_dotplot.pdf"), width=11.6, height=4)
DotPlot(x, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2) + RotatedAxis()
dev.off()

#Do violinplot
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_violinplot.pdf"), width=11.6, height=4)
VlnPlot(x, features = gene_list, split.by = "DEGcluster", stack = TRUE)
dev.off()

#Do heatmap
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_heatmap.pdf"), width=11.6, height=4)
#DefaultAssay(x)<-"RNA"
DoHeatmap(x, features = gene_list, size = 3)
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_heatmaptest.pdf"), width=11.6, height=4)
#DefaultAssay(x)<-"RNA"
DoHeatmap(x, features = gene_list, size = 3)
dev.off()



#ROBO list
list_name<-"ROBO"
gene_list<-read.csv(file = str_glue("{list_dir}/ROBO_intersection_list.csv"))
gene_list<-gene_list$x
gene_list<-na.omit(gene_list)

#Create dotplot
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_dotplot.pdf"), width=11.6, height=4)
DotPlot(x, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2) + RotatedAxis()
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_heatmap.pdf"), width=11.6, height=4)
DoHeatmap(x, features = gene_list, size = 3)
dev.off()

#SLIT list
list_name<-"SLIT"
gene_list<-read.csv(file = str_glue("{list_dir}/SLIT_intersection_list.csv"))
gene_list<-gene_list$x
gene_list<-na.omit(gene_list)

#Create dotplot
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_dotplot.pdf"), width=11.6, height=4)
DotPlot(x, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2) + RotatedAxis()
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_heatmap.pdf"), width=11.6, height=4)
DoHeatmap(x, features = gene_list, size = 3)
dev.off()

#40S
list_name<-"40S"
gene_list<-read.csv(file = str_glue("{list_dir}/40S_intersection_list.csv"))
gene_list<-gene_list$x
gene_list<-na.omit(gene_list)

#Create dotplot
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_dotplot.pdf"), width=11.6, height=4)
DotPlot(x, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2) + RotatedAxis()
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_heatmap.pdf"), width=11.6, height=4)
DoHeatmap(x, features = gene_list, size = 3)
dev.off()

#Peptide
list_name<-"Peptide"
gene_list<-read.csv(file = str_glue("{list_dir}/Peptide_intersection_list.csv"))
gene_list<-gene_list$x
gene_list<-na.omit(gene_list)

#Create dotplot
pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_dotplot.pdf"), width=11.6, height=4)
DotPlot(x, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2) + RotatedAxis()
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_{list_name}_heatmap.pdf"), width=11.6, height=4)
DoHeatmap(x, features = gene_list, size = 3)
dev.off()

#Violin Plot for the 3 differentially expressed genes in all patients
DefaultAssay(x)<-"SCT"
pdf(str_glue("{fig_dir}/{sample_name}_violinplot.pdf"), width=11.6, height=4)
gene_list<- c("RPL5" , "RPS24", "RPL26")
VlnPlot(x, features = gene_list, split.by = "DEGcluster")
dev.off()
#Heatmap for these same three genes
pdf(str_glue("{fig_dir}/{sample_name}_heatmap_3genes.pdf"), width=11.6, height=4)
DoHeatmap(x, features = gene_list, size = 3)
dev.off()
#Do dotplot for the three common genes
pdf(str_glue("{fig_dir}/{sample_name}_dotplot_3genes.pdf"), width=11.6, height=4)
DotPlot(x, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2) + RotatedAxis()
dev.off()

#Do a violin plot for HPV counts by group
pdf(str_glue("{fig_dir}/{sample_name}_violinplot_HPVcounts.pdf"), width=11.6, height=6)
VlnPlot(x, features = "hpv_all", split.by = "DEGcluster", pt.size = 0) + ylim(0, 100) + geom_boxplot(width = 0.1)
dev.off()

#Look at the intersection of the genes in each list because they are quite similar

a<-read.csv(file = str_glue("{list_dir}/viralmRNA_intersection_list.csv"))
a<-a$x
b<-read.csv(file = str_glue("{list_dir}/ROBO_intersection_list.csv"))
b<-b$x
c<-read.csv(file = str_glue("{list_dir}/SLIT_intersection_list.csv"))
c<-c$x
d<-read.csv(file = str_glue("{list_dir}/40S_intersection_list.csv"))
d<-d$x
e<-read.csv(file = str_glue("{list_dir}/Peptide_intersection_list.csv"))
e<-e$x

#Find the intersection between all the lists
peplist<-Reduce(intersect, list(a,b,c,d,e))
#"RPL5"  "RPS24" "RPL26"
write.csv(peplist, file = str_glue("{list_dir}/common_genes_across_pathways.csv"))

robolist<-Reduce(intersect, list (b,c))


#SOOO..I have curated a list from all the enriched pathways:

gene_list<-c("RHOA", "SLIT1", "ROBO2", "EIF2A", "EIF4E", "RPS10", "RPL26", "RPS6", "RPL5", "RPL11", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10")
pdf(str_glue("{fig_dir}/{sample_name}_dotplot_curatedlist.pdf"), width=4, height=8.2)
DotPlot_scCustom(x, features = gene_list, scale.min = 0.8, scale.max = 1, col.min = -2, col.max =2, flip_axes = T, colors_use = c("grey", "purple")) + RotatedAxis()
dev.off()
