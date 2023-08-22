#I have already applied azimuth to the immune cells and the objects are here:
#/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data
#I now need to tweak the annotation a bit
#Need to put in activated T cells
#Need to put in macrophages
#Pull out vectors with cell IDs and immune cell annotation to put into main objects

#fpga
#until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
#conda activate r_4
#cd /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data


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
output<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data"

samples <- c("LN14_filter_SCT_PCA_simpleannotation", "LN23_rep1_filter_SCT_PCA_simpleannotation", "LN23_rep2_filter_SCT_PCA_simpleannotation", "pri34_filter_SCT_PCA_simpleannotation", "pri234_filter_SCT_PCA_simpleannotation", "RZ726_LN2_filter_SCT_PCA_simpleannotation", "RZLN_filter_SCT_PCA_simpleannotation","RZpri_filter_SCT_PCA_simpleannotation")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("read in object")
immune<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}immune_azimuth_annotation.rds'))


#celltypes<-(table(immune$predicted.celltype.l2))
#print("write cell types table")
#write.csv(celltypes, file = (str_glue("{data_dir}/{sample_name}/Initial_cell_types.csv")))

#This reads in the csv file without any of the extra crap
#celltypes<-read.csv(file = "Initial_cell_types.csv")[ ,2:3]

print("create monocyte object")
#I think in order to find macrophages I'll try this:
mono<-subset(immune, subset = predicted.celltype.l2 == c("CD14 Mono", "cDC1", "cDC2", "pDC"))
#works.

print("create counts to assign macrophages")
#I don't think I'll need to get rid of any of the cell types.  Having gone through them I'll need to add in macrophages though
#This is previous code I have used:
#Now try and work out the monocyte population
#Firstly use FCGR2A for macrophages
FCGR2A <-  c("FCGR2A")
MacroCounts <- mono@assays$RNA@counts[which(row.names(mono@assays$RNA@counts) %in%FCGR2A),] %>% as.matrix
# Add HPV expression as a metadata column
mono$MacroCounts <- MacroCounts

print("assign the mean")
meanmacro<-mean(mono$MacroCounts)
print("assign the standard deviation")
sdmacro<-sd(mono$MacroCounts)
print("calculate 2xsd above the mean")
bench<-(meanmacro + (2*sdmacro))

print("put macro counts into immune object")
MacroCounts <- immune@assays$RNA@counts[which(row.names(immune@assays$RNA@counts) %in%FCGR2A),] %>% as.matrix
# Add HPV expression as a metadata column
immune$MacroCounts <- MacroCounts
#Now add in the macrophages to the hn object (taking 2 SD above the mean)
immune$final_immune <- ifelse((immune$MacroCounts>bench) & ((as.character(immune$predicted.celltype.l2)=="CD14 Mono") | (as.character(immune$predicted.celltype.l2)=="cDC1") | (as.character(immune$predicted.celltype.l2)=="cDC2") | (as.character(immune$predicted.celltype.l2)=="pDC")), "Macrophage",as.character(immune$predicted.celltype.l2))

print("UMAP printing")
Idents(immune)<-"final_immune"
pdf(str_glue("{fig_dir}/{sample_name}/UMAP_macrophages.pdf"), width=11.6, height=8.2)
DimPlot(immune, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()

print("save out")
saveRDS(immune, file = str_glue("{output}/{sample_name}/{sample_name}_macrophage_immune_azimuth_annotation.rds"))

#Try to find activated t-cells
#First re-cluster the data
immune<-FindNeighbors(immune, dims = 1:10)
immune<-FindClusters(immune, resolution = 0.5)

#UMAP plot
immune<-RunUMAP(immune, dims = 1:10)

#Plot
Idents(immune)<-"final_immune"
pdf("Reclustering_UMAP.pdf", width=11.6, height=8.2)
DimPlot(immune, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()

#Plot

p1<- DimPlot(immune, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(immune, reduction = "umap", group.by = "final_immune", label = TRUE, repel = TRUE) + NoLegend()
pdf("newclustering_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()
#Worked.

#Try some nebulosa plots for the activated t-cells.  Here is waht Jose told me to do:
#plot_density(tcells, "CCR7")

#Hi Venessa, I would look for expression of IFIT1, IFIT3, and ISG15 for interferon induced T cells
#In case you have any
#and CD137 for activated T cells
pdf("InducedTCells_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(immune, c("IFIT1", "IFIT3", "ISG15"))
dev.off()

#Activated T-cells
pdf("Activated_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(immune, "CD137")
dev.off()

pdf("Activated_IFNG_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(immune, "IFNG")
dev.off()

pdf("Activated_HAVCR2_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(immune, "HAVCR2")
dev.off()

#NK cells
pdf("NK_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(immune, c("CD3D", "NCAM1"))
dev.off()

#Maybe I need to subset out the t cells to see if I can get further granularity
tcell<-subset(immune, subset = ((final_immune == "CD4 Naive") | (final_immune == "CD8 Naive") | (final_immune == "Treg") | (final_immune == "CD4 Proliferating") | (final_immune == "CD8 Proliferating") | (final_immune == "CD4 TCM") | (final_immune == "CD8 TCM") | (final_immune =="dnT") | (final_immune == "CD4 TEM") | (final_immune == "CD8 TEM")))

#Recluster and map
tcell<-FindNeighbors(tcell, dims = 1:10)
tcell<-FindClusters(tcell, resolution = 0.4)

#UMAP plot
tcell<-RunUMAP(tcell, dims = 1:10)

#Plot
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(tcell, reduction = "umap", group.by = "final_immune", label = TRUE, repel = TRUE) + NoLegend()
pdf("newclustering_t_cell_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Try t cell plots again
pdf("InducedTCells_tcell_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(tcell, c("IFIT1", "IFIT3", "ISG15"))
dev.off()

#Activated T-cells
pdf("Activated_UMAP_tcell.pdf", width=11.6, height=8.2)
FeaturePlot(tcell, "CD137")
dev.off()

pdf("Activated_IFNG_tcell_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(tcell, "IFNG")
dev.off()

pdf("Activated_HAVCR2_tcell_UMAP.pdf", width=11.6, height=8.2)
FeaturePlot(tcell, "HAVCR2")
dev.off()

#Looks like I can tell the cytotoxic t cells






