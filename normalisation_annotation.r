#Need to take all the pooled objects normalise them and generate plots so I can roughly annotate them.

# Load required R packages
suppressPackageStartupMessages({
library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(glue)
})

#Define directories
data_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data"
fig_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures"
code_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"


samples <- c("pri234", "pri34", "RZLN", "RZpri", "LN14", "LN23_rep1", "LN23_rep2", "RZ726_LN2")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

#Create figure directories:
if(!dir.exists(str_glue('{fig_dir}/{sample_name}'))){
  dir.create(str_glue('{fig_dir}/{sample_name}'))
}

#Load data
print("load data")
HN<-readRDS(str_glue('{data_dir}/{sample_name}_freemux_MT_filtering_viralreads.rds'))

#Normalise (Check this step with Jose and Walter before running - new SCtransform method on the seurat vignette
print("SCtransform")
HN <- SCTransform(HN, method = "glmGamPoi", verbose = FALSE)
HN <- RunPCA(object = HN, npcs = 100, verbose = FALSE)

print("Elbow plot")
pdf(str_glue("{fig_dir}/{sample_name}/elbowplot.pdf"), width=11.6, height=8.2)
ElbowPlot(HN)
dev.off()

#Run PCA and UMAP
print("Run PCA and UMAP")
HN<-RunUMAP(HN, reduction = "pca", dims = 1:15)
HN<-FindNeighbors(HN, reduction= "pca", dims = 1:15)
HN <- FindClusters(HN, resolution = 0.5)



#Save
print("save out")
saveRDS(HN, file = str_glue("{data_dir}/{sample_name}_filter_SCT_PCA.rds"))

#Plot the cell clusters
print("Dimplots for unlabeled UMAPs")
p1<- DimPlot(HN, reduction = "umap", group.by="MajoritySinglet_Individual_Assignment")
p2<-DimPlot(HN, reduction = "umap", label = TRUE)
pdf(str_glue("{fig_dir}/{sample_name}/unlabeledUMAP.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Do basic cell clustering
print("basic clustering plots")
DefaultAssay(HN) <-"SCT"

print("Tcells")
pdf(str_glue("{fig_dir}/{sample_name}/TCells.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("CD3E", "CD4", "CD8A"))
dev.off()

print("Stroma")
pdf(str_glue("{fig_dir}/{sample_name}/stroma.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("COL1A2"))
dev.off()

print("cancer cells")
pdf(str_glue("{fig_dir}/{sample_name}/cancercells.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("EPCAM", "CASP3"))
dev.off()

print("bcells")
pdf(str_glue("{fig_dir}/{sample_name}/bcells.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("MS4A1", "CD79B"))
dev.off()

print("NK")
pdf(str_glue("{fig_dir}/{sample_name}/NK.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("KLRC1", "NCAM1"))
dev.off()

print("monocytes")
pdf(str_glue("{fig_dir}/{sample_name}/monocytes.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("CD14", "FCGR3A"))
dev.off()

print("plasma")
pdf(str_glue("{fig_dir}/{sample_name}/plasma.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("JCHAIN"))
dev.off()

print("DC")
pdf(str_glue("{fig_dir}/{sample_name}/DC.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("SERPINF1"))
dev.off()

print("blood")
pdf(str_glue("{fig_dir}/{sample_name}/blood.pdf"), width=11.6, height=8.2)
FeaturePlot(HN, features =c("VWF"))
dev.off()


