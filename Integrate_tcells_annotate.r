#I think it would be better to integrate the t-cells and then try and work out which ones are the cytotoxic ones.
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

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/figures/Integrated_tcells"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/code"
output<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data"

print("read in t cell objects")

a<-readRDS(str_glue('{data_dir}/LN14_filter_SCT_PCA_simpleannotation/LN14_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
b<-readRDS(str_glue('{data_dir}/LN23_rep1_filter_SCT_PCA_simpleannotation/LN23_rep1_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
c<-readRDS(str_glue('{data_dir}/LN23_rep2_filter_SCT_PCA_simpleannotation/LN23_rep2_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
d<-readRDS(str_glue('{data_dir}/pri34_filter_SCT_PCA_simpleannotation/pri34_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
e<-readRDS(str_glue('{data_dir}/pri234_filter_SCT_PCA_simpleannotation/pri234_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
f<-readRDS(str_glue('{data_dir}/RZ726_LN2_filter_SCT_PCA_simpleannotation/RZ726_LN2_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
g<-readRDS(str_glue('{data_dir}/RZLN_filter_SCT_PCA_simpleannotation/RZLN_filter_SCT_PCA_simpleannotation_tcell_only.rds'))
h<-readRDS(str_glue('{data_dir}/RZpri_filter_SCT_PCA_simpleannotation/RZpri_filter_SCT_PCA_simpleannotation_tcell_only.rds'))

print("combine")
comb<- list(a, b, c, d, e, f, g, h)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b, c, d, e, f, g, h), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b, c, d, e, f, g, h), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

print("PCA and UMAP")
#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.4)

print("plot")
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(exp.integrated, reduction = "umap", group.by = "final_immune", label = TRUE, repel = TRUE) + NoLegend()
pdf(str_glue("{fig_dir}/UMAP_tcells_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

print("plot by individual")
pdf(str_glue("{fig_dir}/UMAP_by_individual.pdf"), width=11.6, height=8.2)
DimPlot(exp.integrated, reduction = "umap", group.by="seurat_clusters", split.by = "MajoritySinglet_Individual_Assignment", label = TRUE, repel = TRUE) + NoLegend()
dev.off()


print("Print featureplots")
print("HAVCR2 and IFNG")
pdf(str_glue("{fig_dir}/Featureplot_HAVCR2_IFNG.pdf"), width=11.6, height=8.2)
FeaturePlot(exp.integrated, c("HAVCR2", "IFNG"))
dev.off()

print("NK")
#NK cells
pdf(str_glue("{fig_dir}/Featureplot_NK.pdf"), width=11.6, height=8.2)
FeaturePlot(exp.integrated, c("CD3D", "NCAM1"))
dev.off()

print("induced t cells")
pdf(str_glue("{fig_dir}/Featureplot_InducedTcell.pdf"), width=11.6, height=8.2)
FeaturePlot(exp.integrated, c("IFIT1", "IFIT3", "ISG15"))
dev.off()

print("CD137")
#Activated T-cells
pdf(str_glue("{fig_dir}/Featureplot_CD137.pdf"), width=11.6, height=8.2)
FeaturePlot(exp.integrated, "CD137")
dev.off()

print("IFNG")
pdf(str_glue("{fig_dir}/Featureplot_IFNG.pdf"), width=11.6, height=8.2)
FeaturePlot(exp.integrated, "IFNG")
dev.off()

print("HAVCR2")
pdf(str_glue("{fig_dir}/Featureplot_HAVCR2.pdf"), width=11.6, height=8.2)
FeaturePlot(exp.integrated, "HAVCR2")
dev.off()

print("saveout")
saveRDS(exp.integrated, file = str_glue("{output}/integrated_tcell_only.rds"))

print("end")
