#Put in previous transcription and clonal groups from the individual objects into the combined object for dotplots
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


#Data is here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data

sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"


#Read object

comb<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/combined_cancerclones_only.rds")
#This is very complex...the object has many different cell ids which is going to make this tricky
Idents(comb)<-"Freemuxlet_Individual_Assignment"
test<-subset(comb, idents = "HN200519A")
test<-subset(comb, idents = "HN120819A")
test<-subset(comb, idents = "HN021219A")
test<-subset(comb, idents = "HN230620A")
#For HN200519A = the prefix is 20_
#For HN120819A the prefix is 12_
#For HN021219A the prefix is 02_
#For HN230620A the prefix is 23_

#Start with HN200519A
sample_name<-"HN200519A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x<-subset(x, idents = c("Lymph_expanding", "Primary_4", "Primary_6"))
#Collapse into Primary and LN
x$combined<-ifelse(((x$DEGcluster == "Primary_4") | (x$DEGcluster == "Primary_6")), "Primary", "LN")
#Extract the cell ids
new<-x$combined
#Create data frame
df<-data.frame(new)
#Add the prefix
rownames(df) = paste0("20_", rownames(df))
pt1<-df

#Now try and insert
comb<-AddMetaData(object = comb, metadata = df)

#Worked
#HN120819A
sample_name<-"HN120819A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x$combined<-ifelse(((x$DEGcluster =="Lymph_2") | (x$DEGcluster =="Lymph_3")), "LN", as.character(x$DEGcluster))
x$combined<-ifelse(((x$DEGcluster =="primary_2") | (x$DEGcluster =="primary_3")), "Primary", as.character(x$combined))
Idents(x)<-"combined"
x<-subset(x, idents = c("Primary", "LN"))
#Extract the cell ids
new<-x$combined
#Create data frame
df<-data.frame(new)
#In the combined object, there is  "20" at the beginning of each cell ID.  I will need to 
rownames(df) = paste0("12_", rownames(df))
pt2<-df

#Now try and insert
comb<-AddMetaData(object = comb, metadata = df)

#HN021219A
sample_name<-"HN021219A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x$combined<-ifelse(((x$DEGcluster =="Lymph_3") | (x$DEGcluster =="Lymph_4")), "LN", as.character(x$DEGcluster))
x$combined<-ifelse(((x$DEGcluster =="primary")), "Primary", as.character(x$combined))
Idents(x)<-"combined"
x<-subset(x, idents = c("Primary", "LN"))
#Extract the cell ids
new<-x$combined
#Create data frame
df<-data.frame(new)
#In the combined object, there is  "20" at the beginning of each cell ID.  I will need to 
rownames(df) = paste0("02_", rownames(df))
pt3<-df

#Now try and insert
comb<-AddMetaData(object = comb, metadata = df)

#HN230620A
sample_name<-"HN230620A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x$combined<-ifelse(((x$DEGcluster =="Lymph_node")), "LN", as.character(x$DEGcluster))
x$combined<-ifelse(((x$DEGcluster =="primary")), "Primary", as.character(x$combined))
Idents(x)<-"combined"
x<-subset(x, idents = c("Primary", "LN"))
#Extract the cell ids
new<-x$combined
#Create data frame
df<-data.frame(new)
#In the combined object, there is  "20" at the beginning of each cell ID.  I will need to 
rownames(df) = paste0("23_", rownames(df))
pt4<-df

#Try and rbind the tables together
test<-do.call(rbind, list(pt1, pt2, pt3, pt4))

#Now try and insert
comb<-AddMetaData(object = comb, metadata = test)

#No variable features
comb <- SCTransform(comb, vars.to.regress = "percent.mt", verbose = FALSE)
comb <- FindVariableFeatures(comb, selection.method = "vst", nfeatures = 2000)
comb <- RunPCA(comb, features = VariableFeatures(object = comb))

#Save out
saveRDS(comb, file = str_glue("{data_dir}/combined_cancer_only_DEG_transcriptional_groups.rds"))

#Now after that ordeal, need to do a dotplot

#Subset the matrix to take out the NA values
Idents(comb)<-"new"
test<-subset(comb, idents = c("Primary", "LN"))

#Use "rna" slot
DefaultAssay(test)<-"SCT"
gene_list<-c("RHOA", "SLIT1", "ROBO2", "EIF2A", "EIF4E", "RPS10", "RPL26", "RPS6", "RPL5", "RPL11", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10")
pdf(str_glue("{fig_dir}/combined_dotplot_curatedlist_sct.pdf"), width=4, height=8.2)
DotPlot_scCustom(test, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis()
dev.off()

DefaultAssay(test)<-"RNA"
gene_list<-c("RHOA"  , "SLIT1",  "ROBO2" , "EIF2A"  ,"EIF4E",  "RPS10" , "RPL26" , "RPS6" , 
"RPL5" ,  "RPL11" , "PSMB5" , "PSMB6" , "PSMB7" , "PSMB8" , "PSMB9" , "PSMB10",
"RPL10" , "RPL12" , "RPL21" , "RPS24" , "RPL3" ,  "RPL30",  "RPL6" ,  "RPL7"  ,
"RPLP1" , "RPS3" ,  "RPS4X" , "RPS8"  , "EIF3A" , "EIF3E" , "ROCK1" )
pdf(str_glue("{fig_dir}/combined_dotplot_EXTENDEDLIST_sct.pdf"), width=4, height=8.2)
DotPlot_scCustom(test, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis()
dev.off()

# I have a brand new combined object (see other code)


#Use "rna" slot

#Reorder the identities:
Idents(merged)<-"final"
merged@active.ident <- factor(merged@active.ident, 
                            levels=c("Primary", "LN"))

DefaultAssay(merged)<-"SCT"
gene_list<-c("RHOA", "SLIT1", "ROBO2", "EIF2A", "EIF4E", "RPS10", "RPL26", "RPS6", "RPL5", "RPL11", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10")
pdf(str_glue("{fig_dir}/newcombined_object_dotplot_sct.pdf"), width=4, height=8.2)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -2, col.max =2, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis()
dev.off()

DefaultAssay(merged)<-"RNA"
pdf(str_glue("{fig_dir}/newcombined_object_dotplot_rna.pdf"), width=4, height=8.2)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis()
dev.off()

DefaultAssay(merged)<-"SCT"
gene_list<-c("RHOA"  , "SLIT1",  "ROBO2" , "EIF2A"  ,"EIF4E",  "RPS10" , "RPL26" , "RPS6" , 
"RPL5" ,  "RPL11" , "PSMB5" , "PSMB6" , "PSMB7" , "PSMB8" , "PSMB9" , "PSMB10",
"RPL10" , "RPL12" , "RPL21" , "RPS24" , "RPL3" ,  "RPL30",  "RPL6" ,  "RPL7"  ,
"RPLP1" , "RPS3" ,  "RPS4X" , "RPS8"  , "EIF3A" , "EIF3E" , "ROCK1" )
pdf(str_glue("{fig_dir}/new_combined_object_EXTENDEDLIST_sct.pdf"), width=4, height=8.2)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis()
dev.off()

DefaultAssay(merged)<-"RNA"
gene_list<-c("RHOA"  , "SLIT1",  "ROBO2" , "EIF2A"  ,"EIF4E",  "RPS10" , "RPL26" , "RPS6" , 
"RPL5" ,  "RPL11" , "PSMB5" , "PSMB6" , "PSMB7" , "PSMB8" , "PSMB9" , "PSMB10",
"RPL10" , "RPL12" , "RPL21" , "RPS24" , "RPL3" ,  "RPL30",  "RPL6" ,  "RPL7"  ,
"RPLP1" , "RPS3" ,  "RPS4X" , "RPS8"  , "EIF3A" , "EIF3E" , "ROCK1" )
pdf(str_glue("{fig_dir}/new_combined_object_EXTENDEDLIST_rna.pdf"), width=4, height=8.2)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0.5, scale.max = 1, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis()
dev.off()

#Try a heat map
DefaultAssay(merged)<-"SCT"
gene_list<-c("RHOA"  , "SLIT1",  "ROBO2" , "EIF2A"  ,"EIF4E",  "RPS10" , "RPL26" , "RPS6" , 
"RPL5" ,  "RPL11" , "PSMB5" , "PSMB6" , "PSMB7" , "PSMB8" , "PSMB9" , "PSMB10",
"RPL10" , "RPL12" , "RPL21" , "RPS24" , "RPL3" ,  "RPL30",  "RPL6" ,  "RPL7"  ,
"RPLP1" , "RPS3" ,  "RPS4X" , "RPS8"  , "EIF3A" , "EIF3E" , "ROCK1" )
pdf(str_glue("{fig_dir}/new_combined_object_heatmap.pdf"), width=4, height=8.2)
DoHeatmap(merged, features = gene_list) +NoLegend()
dev.off()


#Ii have done more research and now have a bit more of an idea of what I'm looking for.
#Use the combined object for now
#/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/merged_object_paper_SCT.rds
merged<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/merged_object_paper_SCT.rds")

#Calculate the percentage ribosomal genes
#this is wrong:  merged$percent.rb<-PercentageFeatureSet(merged, pattern = "^RB")
merged$percent.ribo<-PercentageFeatureSet(merged, pattern = "^RP[SL]")

#Do a violin plot which is split by primary and LN

Idents(merged)<-"final"
merged@active.ident <- factor(merged@active.ident, 
                            levels=c("Primary", "LN"))
pdf(str_glue("{fig_dir}/new_combined_vln_ribo_genes.pdf"), width=8.2, height=11.6)
p<-VlnPlot(merged, features = c("percent.ribo"), split.by = "final") 
p  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
geom_signif(comparisons = list(c("Primary", "LN")),
              map_signif_level = function(p) sprintf("*p = %.2g", p),
              step_increase = 0.15)

dev.off()

#Do Violinplot for HLA
Idents(merged)<-"final"
merged@active.ident <- factor(merged@active.ident, 
                            levels=c("Primary", "LN"))
pdf(str_glue("{fig_dir}/new_combined_MHCI.pdf"), width=8.2, height=11.6)
p<-VlnPlot(merged, features = c("HLA-A", "HLA-B", "HLA-C"), split.by = "MajoritySinglet_Individual_Assignment") 
p  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
geom_signif(comparisons = list(c("Primary", "LN")),
              map_signif_level = function(p) sprintf("*p = %.2g", p),
              step_increase = 0.15)

dev.off()

my_comparisons<-list(c("Primary", "LN"))
pdf(str_glue("{fig_dir}/new_combined_vln_ribo_genes_ttest.pdf"), width=8.2, height=20)
p<-VlnPlot(merged, features = c("percent.ribo"), split.by = "final") 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + ylim(0, 50)
  dev.off()


my_comparisons<-list(c("Primary", "LN"))
pdf(str_glue("{fig_dir}/new_combined_vln_ribo_genes_wilcox.pdf"), width=8.2, height=11.6)
p<-VlnPlot(merged, features = c("percent.ribo"), split.by = "final") 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons)
  dev.off()

#Nicer graphic
pdf(str_glue("{fig_dir}/new_combined_vln_ribo_genes_nostats.pdf"), width=6, height=4)
p<-VlnPlot(merged, features = c("percent.ribo"), split.by = "final") 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +stat_compare_means(comparisons = my_comparisons)  + NoLegend()
  dev.off()

#No legends

pdf(str_glue("{fig_dir}/new_combined_vln_ribo_genes_nostats_nolegend.pdf"), width=6, height=4)
p<-VlnPlot(merged, features = c("percent.ribo"), split.by = "final", pt.size = 0.1) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons)  + NoLegend() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + ylim(0, 45)+
stat_summary(fun = median, geom='point', size = 40, colour = "black", shape = 95)
  dev.off()


#Do the same for HPV counts

pdf(str_glue("{fig_dir}/new_combined_vln_hpvcounts_genes_nostats.pdf"), width=6, height=4)
p<-VlnPlot(merged, features = c("hpv_all"), split.by = "final") 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons) + ylim(0,1.8) + NoLegend()
  dev.off()
  #Worked and is statistically different!

  #Now hopefully this dot plot looks good:
DefaultAssay(merged)<-"SCT"
  gene_list<-c("TP53", "BAX", "FAS", "BAK1","UBE3A","EIF4E", "EIF2S1","EIF2AK2","PPP1CC","PPP1CB","PPP1CA","PPP1R15A")

pdf(str_glue("{fig_dir}/new_combined_object_Protein_production_dotplot_nolegends.pdf"), width=1.5, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()


pdf(str_glue("{fig_dir}/new_combined_object_Protein_production_dotplot_legends.pdf"), width=4, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, col.min = -4, col.max =4,flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()


#Now do an immune plot
DefaultAssay(merged)<-"SCT"
gene_list<-c( "STAT2", "STAT1", "JAK2", "JAK1", "TYK2",  "IFNAR2","IFNAR1", "IFNGR2","IFNGR1", "IFNLR1" ,"IRF3",  "IRF7","NFKB1","IFIH1")

pdf(str_glue("{fig_dir}/new_combined_object_immune_evation_initial_nolegends.pdf"), width=1.5, height=6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_immune_evation_initial_legends.pdf"), width=4, height=6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

#Downregulation of EIF2alpha (EIF2S1) and PKR (P)

#MHC class III https://www.wikidoc.org/index.php/Major_histocompatibility_complex_class_III_gene_family
#ISG transcription https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7127685/

DefaultAssay(merged)<-"SCT"
gene_list<-c("TNF", "LTA", "LTB", "BAK1", "CFB", "STK19", "HSPA1A", "HSP90AA1", "HSP90B1", "HLA-C", "HLA-B", "HLA-A",  "MX2", "IFITM3", "ISG15","OASL")
pdf(str_glue("{fig_dir}/new_combined_object_immune_evation_downstream_nolegends.pdf"), width=1.5, height=7)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

gene_list<-c("TNF", "LTA", "LTB", "BAK1", "CFB", "STK19", "HSPA1A", "HSP90AA1", "HSP90B1", "HLA-C", "HLA-B", "HLA-A",  "MX2", "IFITM3", "ISG15","OASL")
pdf(str_glue("{fig_dir}/new_combined_object_immune_evation_downstream_legends.pdf"), width=4, height=8)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, col.min = -4, col.max =4, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

# this is a key paper:  https://www.cell.com/cell-reports/pdf/S2211-1247(23)00519-3.pdf

#I think we need a protein list (done), type I IFN inducible genes, then antigen processing pathways.

#IFN inducible genes
DefaultAssay(merged)<-"SCT"
gene_list<-c("IFIT1", "IFIT3", "ISG15", "MX1", "MX2", "IRF3", "IRF7", "JAK1", "JAK2", "STAT1", "STAT2", "TYK2")
pdf(str_glue("{fig_dir}/new_combined_object_IFN_INDUCIBLE_nolegends.pdf"), width=1.5, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

gene_list<-c("IFIT1", "IFIT3", "ISG15", "MX1", "MX2", "IRF3", "IRF7", "JAK1", "JAK2", "STAT1", "STAT2", "TYK2")
pdf(str_glue("{fig_dir}/new_combined_object_IFN_INDUCIBLE_legends.pdf"), width=4, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()


#Antigen processing pathways
DefaultAssay(merged)<-"SCT"
gene_list<-c("CFB","HSP90B1","HSPA1A", "HLA-DQB1","HLA-DQA1","HLA-DRB1", "HLA-DRA", "HLA-C", "HLA-B", "HLA-A", "B2M", "TAPBP", "TAP2", "TAP1", "PSMB9", "PSMB8", "STING","IFIH1")
pdf(str_glue("{fig_dir}/new_combined_object_IMMUNE_PTH_nolegends.pdf"), width=1.5, height=6.5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

gene_list<-c("CFB","HSP90B1","HSPA1A", "HLA-DQB1","HLA-DQA1","HLA-DRB1", "HLA-DRA", "HLA-C", "HLA-B", "HLA-A", "B2M", "TAPBP", "TAP2", "TAP1", "PSMB9", "PSMB8", "STING","IFIH1")
pdf(str_glue("{fig_dir}/new_combined_object_IMMUNE_PTH_legends.pdf"), width=4, height=6.5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

#Try another round because increased experession of MHCI molecules doesn't make a lot of sends
DefaultAssay(merged)<-"SCT"
gene_list<-c("HLA-C", "HLA-B", "HLA-A", "B2M", "TAP1", "TAP2", "TAPBP", "HLA-DRB1", "HLA-DRA", "DLA-DRB3", "HLA-DRB4", "HLA-DQB1", "HLA-DQA1")
pdf(str_glue("{fig_dir}/new_combined_object_MHCII_nolegends.pdf"), width=1.5, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

gene_list<-c("HLA-C", "HLA-B", "HLA-A", "B2M", "TAP1", "TAP2", "TAPBP", "HLA-DRB1", "HLA-DRA", "DLA-DRB3", "HLA-DRB4", "HLA-DQB1", "HLA-DQA1")
pdf(str_glue("{fig_dir}/new_combined_object_MHCII_legends.pdf"), width=4, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

#Try some more translation initiation factors
DefaultAssay(merged)<-"SCT"
gene_list<-c("EIF4E", "EIF4G", "EIF4A","EIF2", "EIF3", "EIF2AK2", "EIF2S1")
pdf(str_glue("{fig_dir}/new_combined_object_40S.pdf"), width=1.5, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

gene_list<-c("EIF4E", "EIF4G", "EIF4A","EIF2", "EIF3", "EIF2AK2", "EIF2S1")
pdf(str_glue("{fig_dir}/new_combined_object_40Snolegends.pdf"), width=4, height=5)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

#God help me...seriously.  https://www.annualreviews.org/doi/pdf/10.1146/annurev-cancerbio-030419-033420

DefaultAssay(merged)<-"SCT"
gene_list<-c("EIF4E", "EIF4G", "EIF4A","EIF2", "EIF3", "EIF2AK2", "EIF2S1", "PDCD4", "RPL36A", "RPL23A", "RL35A", "RPL10", "RPL5", "RPL39","RPL22", "RPS15", "RPS17", "RPS24", "RPS19", "MYC", "RB1", "MTOR")
pdf(str_glue("{fig_dir}/new_combined_cancer_translation.pdf"), width=1.5, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_cancertranslation_nolegends.pdf"), width=4, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()


#After speaking with Erin.  There is CAP dependent translation which is the "normal" method of translation.  Then there is CAP independent translation which is often a stress response
#Try these cap dependent regulatory genes

DefaultAssay(merged)<-"SCT"
gene_list<-c("EIF2AK3", "EIF2AK1", "EIF2AK2", "EIF2AK4", "EIF2S1", "PPP1R15A", "PPP1CA", "PPP1CB", "PPP1CC", "MTOR","ATP7A", "EEF2", "RPS6", "PDCD4", "EIF4EBP1", "EEF2K", "MYC", "TP53", "RB1")
pdf(str_glue("{fig_dir}/new_combined_cancer_translationa_important.pdf"), width=1.5, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_translationa_important_nolegends.pdf"), width=4, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()


#Another set of bloody genes. https://link.springer.com/article/10.1007/s00018-016-2428-2.  Now looking at cap-indpendent translation


DefaultAssay(merged)<-"SCT"
gene_list<-c("VEGFA", "HIF1A", "FGF2", "BCL2", "TRPV6", "SLC38A2", "XIAP", "SREBF1", "BIRCR", "GULP", "PCBP1", "PCBP2", "EIF4G1", "EIF4G2", "EIF4G3", "EIF4A1", "EIF4A2", "EIF4A3")
pdf(str_glue("{fig_dir}/new_combined_cancer_cap_indpendent.pdf"), width=1.5, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_cap_indpendent_nolegends.pdf"), width=4, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

#https://www.cell.com/cell-systems/pdf/S2405-4712(21)00252-0.pdf  is a very interesting paper about cap-independent translation
#Try a scatter plot

#That doesn't work very well.  Try something different

#Calculate the percentage of EIF4G1
merged$percent.EIF4G1<-PercentageFeatureSet(merged, pattern = "EIF4G1")
merged$percent.EIF4E<-PercentageFeatureSet(merged, pattern = "EIF4E")
merged$hpv_all<-PercentageFeatureSet(merged, pattern = "^HpV16")
#Subset where these values =0
submerged<-subset(merged, subset = percent.EIF4E != 0)
submerged<-subset(merged, subset = percent.EIF4G1 != 0)
submerged<-subset(submerged, subset = percent.EIF4G1 != 0)
submerged<-subset(submerged, subset = percent.EIF4A1 !=0)

#Then separate the Primary and LN
DefaultAssay(submerged)<-"SCT"
primary<-subset(merged, subset = final =="Primary")
p<-FeatureScatter(object = primary, feature1 = "percent.EIF4G1", feature2 = "percent.EIF4E", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_capI_scatter_primary.pdf"), width=4, height=4)
p
dev.off()

LN<-subset(merged, subset = final =="LN")
p<-FeatureScatter(object = LN, feature1 = "percent.EIF4G1", feature2 = "percent.EIF4E", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_capI_scatter_LN.pdf"), width=4, height=4)
p
dev.off()

primary<-subset(merged, subset = final =="Primary")
p<-FeatureScatter(object = primary, feature1 = "percent.EIF4G1", feature2 = "percent.EIF4A1", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_capI_scatter_primaryA1.pdf"), width=4, height=4)
p
dev.off()

LN<-subset(merged, subset = final =="LN")
p<-FeatureScatter(object = LN, feature1 = "percent.EIF4G1", feature2 = "percent.EIF4A1", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_capI_scatter_LNA1.pdf"), width=4, height=4)
p
dev.off()


p<-FeatureScatter(object = submerged, feature1 = "EIF4G1", feature2 = "EIF4E", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_E.pdf"), width=4, height=4)
p
dev.off()

DefaultAssay(primary)<-"SCT"
p<-FeatureScatter(object = primary, feature1 = "EIF4G1", feature2 = "EIF4E", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_E_primary.pdf"), width=4, height=4)
p
dev.off()

DefaultAssay(LN)<-"SCT"
p<-FeatureScatter(object = LN, feature1 = "EIF4G1", feature2 = "EIF4E", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_E_LN.pdf"), width=4, height=4)
p
dev.off()

p<-FeatureScatter(object = primary, feature1 = "EIF4G1", feature2 = "EIF4A1", slot = "data", jitter = TRUE)
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_A1_primary.pdf"), width=4, height=4)
p
dev.off()

p<-FeatureScatter(object = LN, feature1 = "EIF4G1", feature2 = "EIF4A1", slot = "data", jitter = TRUE)
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_A1_LN.pdf"), width=4, height=4)
p
dev.off()

p<-FeatureScatter(object = primary, feature1 = "EIF4G1", feature2 = "HpV16gp1", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_E6_primary.pdf"), width=4, height=4)
p
dev.off()

p<-FeatureScatter(object = LN, feature1 = "EIF4G1", feature2 = "HpV16gp1", slot = "data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_v_E6_LN.pdf"), width=4, height=4)
p
dev.off()


#Try with all a different scale methods;
allmerged<-SCTransform(merged, return.only.var.genes = FALSE, vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(allmerged, file = str_glue("{data_dir}/combined_cancer_allmerged_allfeatures.rds"))
primary<-subset(allmerged, subset = final == "Primary")
DefaultAssay(primary)<-"SCT"
p<-FeatureScatter(object = primary, feature1 = "EIF4G1", feature2 = "EIF4E", slot = "scale.data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_4E_primary.pdf"), width=4, height=4)
p
dev.off()

DefaultAssay(primary)<-"SCT"
p<-FeatureScatter(object = primary, feature1 = "EIF4G1", feature2 = "EIF4A1", slot = "scale.data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_A1_primary.pdf"), width=4, height=4)
p
dev.off()

LN<-subset(allmerged, subset = final == "LN")
DefaultAssay(LN)<-"SCT"
p<-FeatureScatter(object = LN, feature1 = "EIF4G1", feature2 = "EIF4E", slot = "scale.data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_4E_LN.pdf"), width=4, height=4)
p
dev.off()


p<-FeatureScatter(object = LN, feature1 = "EIF4G1", feature2 = "EIF4A1", slot = "scale.data")
pdf(str_glue("{fig_dir}/new_combined_object_G1_A1_LN.pdf"), width=4, height=4)
p
dev.off()

DefaultAssay(allmerged)<-"SCT"
gene_list<-c("VEGFA", "HIF1A", "FGF2", "BCL2", "TRPV6", "SLC38A2", "XIAP", "SREBF1", "BIRCR", "GULP", "PCBP1", "PCBP2", "EIF4G1", "EIF4G2", "EIF4G3", "EIF4A1", "EIF4A2", "EIF4A3")
pdf(str_glue("{fig_dir}/new_combined_cancer_cap_indpendent_test.pdf"), width=1.5, height=11.6)
DotPlot_scCustom(allmerged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_cap_indpendent_nolegends_test.pdf"), width=4, height=11.6)
DotPlot_scCustom(allmerged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

#I think I've finally made my mind up

DefaultAssay(merged)<-"SCT"
gene_list<-c("EIF4G1", "EIF4A1", "EIF4E")
pdf(str_glue("{fig_dir}/new_combined_cancer_initiatingcomplex_nolegends.pdf"), width=1.5, height=4)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_initiatingcomplex_legends.pdf"), width=4, height=4)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()




pdf(str_glue("{fig_dir}/new_combined_object_capdependent_legends.pdf"), width=4, height=8)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

gene_list<-c("RPS2", "ERCC6L", "CKAP2", "CCNA2", "MCM7", "VEGFA", "HIF1A", "SP1", "PCBP1", "PCBP2")
pdf(str_glue("{fig_dir}/new_combined_cancer_capINdependent_nolegends.pdf"), width=1.5, height=8)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_capINdependent_legends.pdf"), width=4, height=8)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()

gene_list<-c("HpV16gp1", "HpV16gp2", "HpV16gp3", "HpV16gp4", "HpV16gp5", "HpV16gp6", "HpV16gp7","HpV16gp8")
pdf(str_glue("{fig_dir}/new_combined_cancer_HPV_nolegends.pdf"), width=1.5, height=8)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

pdf(str_glue("{fig_dir}/new_combined_object_HPV_legends.pdf"), width=4, height=8)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue")) + RotatedAxis() 
dev.off()


#Do plots with matri with returning all genes
Idents(allmerged)<-"MajoritySinglet_Individual_Assignment"
gene_list<-c("EIF4G1", "EIF4A1", "EIF4E","EIF4EBP1", "PDCD4", "EIF2S1", "EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "PPP1R15A", "PPP1CA")
pdf(str_glue("{fig_dir}/new_combined_cancer_capdependent_nolegends_allvariablegenes.pdf"), width=1.5, height=8)
DotPlot_scCustom(allmerged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue"), scale = FALSE) + RotatedAxis() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + NoLegend()
dev.off()

gene_list<-c("EIF4G1", "EIF4A1", "EIF4E","EIF4EBP1", "PDCD4", "EIF2S1", "EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "PPP1R15A", "PPP1CA")
pdf(str_glue("{fig_dir}/new_combined_cancer_capdependent_nolegends_allvariablegenes.pdf"), width=4, height=8)
DotPlot_scCustom(allmerged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue"), scale = FALSE) + RotatedAxis()
dev.off()

#Okay, I think the dotplots are fucked.

#This is so depressing.  I'll try a heatmap and some violin plots, but if they don't work, we are going to have to go back to the drawing board.
Idents(allmerged)<-"final"
DefaultAssay(merged)<-"SCT"
gene_list<-c("EIF4E", "EIF4G1", "EIF4A1", "EIF2AK2", "EIF2S1", "PDCD4", "RPL36A", "RPL23A", "RL35A", "RPL10", "RPL5", "RPL39","RPL22", "RPS15", "RPS17", "RPS24", "RPS19", "MYC", "RB1", "MTOR", "EIF4EBP1", "PDCD4", "EIF2S1", "EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "PPP1R15A", "PPP1CA","VEGFA", "HIF1A", "FGF2", "BCL2", "TRPV6", "SLC38A2", "XIAP", "SREBF1", "BIRCR", "GULP", "PCBP1", "PCBP2", "EIF4G1", "EIF4G2", "EIF4G3", "EIF4A1", "EIF4A2", "EIF4A3","CFB","HSP90B1","HSPA1A", "HLA-DQB1","HLA-DQA1","HLA-DRB1", "HLA-DRA", "HLA-C", "HLA-B", "HLA-A", "B2M", "TAPBP", "TAP2", "TAP1", "PSMB9", "PSMB8", "STING","IFIH1")
pdf(str_glue("{fig_dir}/new_combined_cancer_translation.pdf"), width=4, height=11.6)
DotPlot_scCustom(merged, features = gene_list, scale.min = 0, scale.max = 60, flip_axes = T, colors_use = c("grey", "blue"), scale = FALSE) + RotatedAxis() 
dev.off()

#Try a heat map although I have low expectations
pdf(str_glue("{fig_dir}/new_combined_cancer_translation_heatmap.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = mymarkers, size =4, angle =90) + NoLegend()
p 
dev.off()

Idents(allmerged)<-"final"
markers <- FindMarkers(allmerged, ident.1 = "LN", ident.2 = "Primary", verbose = FALSE)
markers<-FindAllMarkers(allmerged, min.pct = 0.1, logfc.threshold = 0.25)
top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
head(markers, n = 100)

markers <- FindMarkers(allmerged, ident.1 = "LN", ident.2 = "Primary", verbose = FALSE)
test<-filter(markers, p_val_adj < 0.05)
test<-arrange(test, p_val_adj)
#top100
p<-head(test,100)
p<-arrange(p, p_val_adj)
list<-rownames(p)
allmerged@active.ident <- factor(allmerged@active.ident, 
                            levels=c("Primary", "LN"))

#I now have this list which is variable genes intersected with all the genes from the enriched pathways
 [1] "EEF1B2"  "EIF3G"   "EIF3I"   "EIF3K"   "PFN1"    "PSMA1"   "PSMA2"  
 [8] "PSMA3"   "PSMA5"   "PSMA7"   "PSMB1"   "PSMB2"   "PSMB3"   "PSMB4"  
[15] "PSMB5"   "PSMB6"   "PSMB7"   "PSMB9"   "PSMC3"   "PSMC4"   "PSMC5"  
[22] "PSMD1"   "PSMD11"  "PSMD13"  "PSMD2"   "PSMD3"   "PSMD4"   "PSMD8"  
[29] "PSMD9"   "PSME1"   "PSME2"   "RAC1"    "RHOA"    "RPL13"   "RPL13A" 
[36] "RPL14"   "RPL18A"  "RPL19"   "RPL21"   "RPL23A"  "RPL26L1" "RPL29"  
[43] "RPL3"    "RPL30"   "RPL37"   "RPL41"   "RPL5"    "RPL6"    "RPL7A"  
[50] "RPL8"    "RPLP0"   "RPLP1"   "RPS11"   "RPS16"   "RPS21"   "RPS24"  
[57] "RPS25"   "RPS26"   "RPS27A"  "RPS27L"  "RPS3"    "RPS5"    "RPS7"   
[64] "RPS9"    "RPSA"    "SEC11C"  "SEC61B"  "SPCS1"   "SRP9"    "UBA52"  
[71] "UBB"   




pdf(str_glue("{fig_dir}/combined_enriched_variable_heatmap.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = a, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

#IFN Genes
[1] "RPS27A" "UBB"    "UBA52"  "EIF4G2" "ARIH1"  "EIF4A2" "FLNA"   "IRF3"  
[9] "ISG15" 

#Adding in some more genes:
 [1] "RPS27A"   "UBB"      "UBA52"    "EIF4G2"   "ARIH1"    "EIF4A2"  
 [7] "FLNA"     "IRF3"     "ISG15"    "TAP1"     "TAP2"     "CFB"     
[13] "HSP90B1"  "HSPA1A"   "HLA-DQB1" "HLA-DQA1" "HLA-DRB1" "HLA-DRA" 
[19] "HLA-C"    "HLA-B"    "HLA-A"    "B2M"      "TAPBP" 

pdf(str_glue("{fig_dir}/combined_IFN_variable_heatmap.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = d, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()


markers1<- presto::wilcoxauc(allmerged, 'final', assay = 'data')
markers1<- top_markers(markers1, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- markers1 %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
   .[!is.na(.)]

   mat<- allmerged[["RNA"]]@data[all_markers, ] %>% as.matrix()
   mat<- t(scale(t(mat)))
   cluster_anno<- allmerged@meta.data$final

   quantile(mat, c(0.1, 0.95))
   Seurat::PurpleAndYellow()
   col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
pdf(str_glue("{fig_dir}/testheatmap.pdf"), width=11.6, height=8.2)
   p<-Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 4),
        column_title_rot = 90,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)
p 
        dev.off()

#How can I get both the patient ID and the location in the same meta data column?  Who the fuck knows.

allmerged$patients<-ifelse(((allmerged$MajoritySinglet_Individual_Assignment == "HN200519A")), "Patient1", "Nil")
allmerged$patients<-ifelse(((allmerged$MajoritySinglet_Individual_Assignment == "HN120819A")), "Patient2", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$MajoritySinglet_Individual_Assignment == "HN021219A")), "Patient3", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$MajoritySinglet_Individual_Assignment == "HN230620A")), "Patient4", as.character(allmerged$patients))

#okay, worked

allmerged$patients<-ifelse(((allmerged$patients == "Patient1" & allmerged$final == "Primary")), "Patient1_Primary", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient1" & allmerged$final == "LN")), "Patient1_LN", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient2" & allmerged$final == "Primary")), "Patient2_Primary", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient2" & allmerged$final == "LN")), "Patient2_LN", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient3" & allmerged$final == "Primary")), "Patient3_Primary", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient3" & allmerged$final == "LN")), "Patient3_LN", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient4" & allmerged$final == "Primary")), "Patient4_Primary", as.character(allmerged$patients))
allmerged$patients<-ifelse(((allmerged$patients == "Patient4" & allmerged$final == "LN")), "Patient4_LN", as.character(allmerged$patients))

Idents(allmerged)<-"patients"
allmerged@active.ident <- factor(allmerged@active.ident, 
                            levels=c("Patient1_Primary",  "Patient2_Primary", "Patient3_Primary",  "Patient4_Primary",  "Patient1_LN", "Patient2_LN","Patient3_LN","Patient4_LN"))
saveRDS(allmerged, file = str_glue("{data_dir}/combined_cancer_allmerged_allfeatures.rds"))


#Now repeat heatmaps
Idents(allmerged)<-"patients"
pdf(str_glue("{fig_dir}/combined_enriched_variable_heatmap.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = a, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

c<-c("RPS27A", "UBB",    "UBA52" , "EIF4G2" ,"ARIH1" , "EIF4A2" ,"FLNA" ,  "IRF3" , "ISG15") 
pdf(str_glue("{fig_dir}/combined_IFN_variable_heatmap.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = c, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

viral_genes<-c("HpV16gp1", "HpV16gp2", "HpV16gp3", "HpV16gp4", "HpV16gp5", "HpV16gp6", "HpV16gp7", "HpV16gp8")
pdf(str_glue("{fig_dir}/combined_viralgenes_heatmap.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = viral_genes, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()


Idents(allmerged)<-"final"
allmerged@active.ident <- factor(allmerged@active.ident, 
                            levels=c("Primary", "LN"))
pdf(str_glue("{fig_dir}/combined_enriched_variable_heatmap_P_LN.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = a, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

c<-c("RPS27A", "UBB",    "UBA52" , "EIF4G2" ,"ARIH1" , "EIF4A2" ,"FLNA" ,  "IRF3" , "ISG15") 
pdf(str_glue("{fig_dir}/combined_IFN_variable_heatmapP_LN.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = c, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

viral_genes<-c("HpV16gp1", "HpV16gp2", "HpV16gp3", "HpV16gp4", "HpV16gp5", "HpV16gp6", "HpV16gp7", "HpV16gp8")
pdf(str_glue("{fig_dir}/combined_viralgenes_heatmapP_LN.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = viral_genes, size =4) + NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()