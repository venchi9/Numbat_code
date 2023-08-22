#Try to find activated t-cells
#Try and find cytotoxic tcells

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
immune<-readRDS(str_glue('{data_dir}/{sample_name}/{sample_name}_macrophage_immune_azimuth_annotation.rds'))

print("create t-cell subgroup")
tcell<-subset(immune, subset = ((final_immune == "CD4 Naive") | (final_immune == "CD8 Naive") | (final_immune == "Treg") | (final_immune == "CD4 Proliferating") | (final_immune == "CD8 Proliferating") | (final_immune == "CD4 TCM") | (final_immune == "CD8 TCM") | (final_immune =="dnT") | (final_immune == "CD4 TEM") | (final_immune == "CD8 TEM")))

print("recluster and map")
tcell<-FindNeighbors(tcell, dims = 1:10)
tcell<-FindClusters(tcell, resolution = 0.4)

#UMAP plot
tcell<-RunUMAP(tcell, dims = 1:10)

print("plot umap")
#Plot
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-DimPlot(tcell, reduction = "umap", group.by = "final_immune", label = TRUE, repel = TRUE) + NoLegend()
pdf(str_glue("{fig_dir}/{sample_name}/UMAP_tcells_only.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

print("plot by individual")
pdf(str_glue("{fig_dir}/{sample_name}/UMAP_by_individual.pdf"), width=11.6, height=8.2)
DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", split.by = "MajoritySinglet_Individual_Assignment", label = TRUE, repel = TRUE) + NoLegend()
dev.off()


print("Print featureplots")
print("HAVCR2 and IFNG")
pdf(str_glue("{fig_dir}/{sample_name}/Featureplot_HAVCR2_IFNG.pdf"), width=11.6, height=8.2)
FeaturePlot(tcell, c("HAVCR2", "IFNG"))
dev.off()

print("NK")
#NK cells
pdf(str_glue("{fig_dir}/{sample_name}/Featureplot_NK.pdf"), width=11.6, height=8.2)
FeaturePlot(tcell, c("CD3D", "NCAM1"))
dev.off()

print("induced t cells")
pdf(str_glue("{fig_dir}/{sample_name}/Featureplot_InducedTcell.pdf"), width=11.6, height=8.2)
FeaturePlot(tcell, c("IFIT1", "IFIT3", "ISG15"))
dev.off()

print("CD137")
#Activated T-cells
pdf(str_glue("{fig_dir}/{sample_name}/Featureplot_CD137.pdf"), width=11.6, height=8.2)
FeaturePlot(tcell, "CD137")
dev.off()

print("IFNG")
pdf(str_glue("{fig_dir}/{sample_name}/Featureplot_IFNG.pdf"), width=11.6, height=8.2)
FeaturePlot(tcell, "IFNG")
dev.off()

print("saveout")
saveRDS(tcell, file = str_glue("{output}/{sample_name}/{sample_name}_tcell_only.rds"))

print("end")

########################
#Do some further markers
#Integrated t cell object is here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data

#Plot the UMAP for each individual
#First create subgroups
#HN021219A
HN021219A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN021219A")
HN070219A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN070219A")
HN120520A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN120520A")
HN120819A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN120819A")
HN170419A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
HN200519A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN200519A")
HN230620A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN230620A")

fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/figures/Integrated_tcells"
#HN021219A
sample_name<-"HN021219A"
print("plot by individual")
p1<- DimPlot(HN021219A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN021219A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN070219A
sample_name<-"HN070219A"
print("plot by individual")
p1<- DimPlot(HN070219A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN070219A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN120520A
sample_name<-"HN120520A"
print("plot by individual")
p1<- DimPlot(HN120520A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN120520A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN120819A
sample_name<-"HN120819A"
print("plot by individual")
p1<- DimPlot(HN120819A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN120819A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN170419A
sample_name<-"HN170419A"
print("plot by individual")
p1<- DimPlot(HN170419A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN170419A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN200519A
sample_name<-"HN200519A"
print("plot by individual")
p1<- DimPlot(HN200519A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN200519A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN230620A
sample_name<-"HN230620A"
print("plot by individual")
p1<- DimPlot(HN230620A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN230620A, "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individual.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Now look at CD29 as a marker
sample_name<-"HN021219A"
print("plot by individual")
p1<- DimPlot(HN021219A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN021219A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN070219A
sample_name<-"HN070219A"
print("plot by individual")
p1<- DimPlot(HN070219A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN070219A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN120520A
sample_name<-"HN120520A"
print("plot by individual")
p1<- DimPlot(HN120520A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN120520A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN120819A
sample_name<-"HN120819A"
print("plot by individual")
p1<- DimPlot(HN120819A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN120819A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN170419A
sample_name<-"HN170419A"
print("plot by individual")
p1<- DimPlot(HN170419A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN170419A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN200519A
sample_name<-"HN200519A"
print("plot by individual")
p1<- DimPlot(HN200519A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN200519A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN230620A
sample_name<-"HN230620A"
print("plot by individual")
p1<- DimPlot(HN230620A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN230620A, "ITGB1")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualITGB1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Now interleukin 2
sample_name<-"HN021219A"
print("plot by individual")
p1<- DimPlot(HN021219A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN021219A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN070219A
sample_name<-"HN070219A"
print("plot by individual")
p1<- DimPlot(HN070219A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN070219A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN120520A
sample_name<-"HN120520A"
print("plot by individual")
p1<- DimPlot(HN120520A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN120520A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN120819A
sample_name<-"HN120819A"
print("plot by individual")
p1<- DimPlot(HN120819A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN120819A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN170419A
sample_name<-"HN170419A"
print("plot by individual")
p1<- DimPlot(HN170419A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN170419A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN200519A
sample_name<-"HN200519A"
print("plot by individual")
p1<- DimPlot(HN200519A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN200519A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#HN230620A
sample_name<-"HN230620A"
print("plot by individual")
p1<- DimPlot(HN230620A, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-FeaturePlot(HN230620A, "IL2")
pdf(str_glue("{fig_dir}_{sample_name}_UMAP_by_individualIL2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#This is all very puzzling.  I think I will try some ridge plots
sample_name<-"t-cells_integrated"
pdf(str_glue("{fig_dir}_{sample_name}_ridgeplot_cytotoxic_tcells.pdf"), width=11.6, height=8.2)
RidgePlot(tcell, features = c("IL2", "IFNG", "ITGB1", "HAVCR2"), ncol =2)
dev.off()

#Violin plot
sample_name<-"t-cells_integrated"
pdf(str_glue("{fig_dir}_{sample_name}_violinplot_cytotoxic_tcells.pdf"), width=11.6, height=8.2)
VlnPlot(tcell, features = c("IL2", "IFNG", "ITGB1", "HAVCR2"), ncol =2)
dev.off()

pdf(str_glue("{fig_dir}_{sample_name}_violinplot_IL2_tcells.pdf"), width=11.6, height=8.2)
VlnPlot(tcell, features = "IL2")
dev.off()

pdf(str_glue("{fig_dir}_{sample_name}_violinplot_IFNG_tcells.pdf"), width=11.6, height=8.2)
VlnPlot(tcell, features = "IFNG")
dev.off()

pdf(str_glue("{fig_dir}_{sample_name}_violinplot_HAVCR2_tcells.pdf"), width=11.6, height=8.2)
VlnPlot(tcell, features = "HAVCR2")
dev.off()

pdf(str_glue("{fig_dir}_{sample_name}_violinplot_ITGB1_tcells.pdf"), width=11.6, height=8.2)
VlnPlot(tcell, features = "ITGB1")
dev.off()

#Do a dotplot
pdf(str_glue("{fig_dir}_{sample_name}_dotplot_cytotoxic_tcells.pdf"), width=11.6, height=8.2)
DotPlot(tcell, features = c("IL2", "IFNG", "ITGB1", "HAVCR2")) + RotatedAxis()
dev.off()

#Heatmap
features <- c("ITGB1", "GNLY", "GZMK", "TNF", "PRF1", "IL2", "PTGS2", "HAVCR2", "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_heatmap_cytotoxic_tcells.pdf"), width=11.6, height=8.2)
DoHeatmap(subset(tcell, downsample = 500), features = features, size =3)
dev.off()

#Maybe try the heatmap with just the effector T cells
CD8<-subset(tcell, subset = final_immune == c("CD8 TEM"))
pdf(str_glue("{fig_dir}_{sample_name}_heatmap_cytotoxic_CD8.pdf"), width=11.6, height=8.2)
DoHeatmap(CD8, features = features, size =3)
dev.off()

#Now try a scatter plot
Idents(tcell)<-"final_immune"
sample_name<-"integrated_tcells"
pdf(str_glue("{fig_dir}_{sample_name}_scatterplot_IFNG_HAVCR2.pdf"), width=11.6, height=8.2)
FeatureScatter(tcell, feature1 = "IFNG", feature2 = "HAVCR2")
dev.off()

sample_name<-"CD8_ALL"
CD8<-subset(tcell, subset = final_immune == c("CD8 TEM", "CD8 TCM", "CD8 Proliferating", "CD8 Naive"))
pdf(str_glue("{fig_dir}_{sample_name}_scatterplot_IFNG_HAVCR2.pdf"), width=11.6, height=8.2)
FeatureScatter(CD8, feature1 = "IFNG", feature2 = "HAVCR2")
dev.off()

#Try and subtype the CD8 TEM
sample_name<-"CD8_TEM"
CD8<-subset(tcell, subset = final_immune == c("CD8 TEM"))
features <- c("ITGB1", "GNLY", "GZMK", "TNF")
pdf(str_glue("{fig_dir}_{sample_name}_violinplot_cytotoxic_CD8TEM.pdf"), width=8.2, height=11.6)
VlnPlot(CD8, features = features, ncol =2)
dev.off()

features <- c("PRF1", "IL2", "HAVCR2", "IFNG")
pdf(str_glue("{fig_dir}_{sample_name}_violinplot_cytotoxic_CD8TEM2.pdf"), width=8.2, height=11.6)
VlnPlot(CD8, features = features, ncol =2)
dev.off()

#Try exhaustion markers
features<-c("PDCD1", "HAVCR2")
pdf(str_glue("{fig_dir}_{sample_name}_violinplot_exhausted_CD8TEM1.pdf"), width=8.2, height=11.6)
VlnPlot(CD8, features = features, ncol =2)
dev.off()


features<-c("LAG3", "BTLA", "CTLA4")
pdf(str_glue("{fig_dir}_{sample_name}_violinplot_exhausted_CD8TEM2.pdf"), width=8.2, height=11.6)
VlnPlot(CD8, features = features, ncol =2)
dev.off()

#Try ridge plots again
sample_name<-"CD8TEM"
pdf(str_glue("{fig_dir}_{sample_name}_ridgeplot_exhaustion_tcells.pdf"), width=11.6, height=8.2)
RidgePlot(CD8, features = c("PDCD1", "HAVCR2","LAG3", "BTLA", "CTLA4"), ncol =2)
dev.off()

sample_name<-"CD8TEM"
pdf(str_glue("{fig_dir}_{sample_name}_ridgeplot_cytotoxic_CD8.pdf"), width=11.6, height=8.2)
RidgePlot(CD8, features = c("GNLY", "GZMK", "TNF"), ncol =2)
dev.off()

sample_name<-"CD8TEM"
pdf(str_glue("{fig_dir}_{sample_name}_ridgeplot_cytotoxic_CD8_2.pdf"), width=11.6, height=8.2)
RidgePlot(CD8, features = c("FASLG", "PRF1", "IL2", "IFNG"), ncol =2)
dev.off()

#Feature scatter with the more representative genes
pdf(str_glue("{fig_dir}_{sample_name}_scatterplot_IFNG_LAG3_CD8TEM.pdf"), width=11.6, height=8.2)
FeatureScatter(CD8, feature1 = "IFNG", feature2 = "LAG3")
dev.off()

pdf(str_glue("{fig_dir}_{sample_name}_scatterplot_PDCD1_LAG3_CD8TEM.pdf"), width=11.6, height=8.2)
FeatureScatter(CD8, feature1 = "PDCD1", feature2 = "LAG3")
dev.off()

#Heatmap with all the exhaustion markers
Idents(CD8)<-"seurat_clusters"
pdf(str_glue("{fig_dir}_{sample_name}_heatmap_exhaustionmarkers_CD8.pdf"), width=11.6, height=8.2)
features<-c("LAG3", "BTLA", "CTLA4","PDCD1", "HAVCR2")
DoHeatmap(CD8, features = features, size =3)
dev.off()

Idents(CD8)<-"seurat_clusters"
pdf(str_glue("{fig_dir}_{sample_name}_heatmap_cytotoxicmarkers_CD8.pdf"), width=11.6, height=8.2)
features<-c("GNLY", "GZMK", "TNF", "FASLG", "PRF1", "IL2", "IFNG")
DoHeatmap(CD8, features = features, size =3)
dev.off()

#These heatmaps are the most helpful.  I will repeat now with all the CD8 cells (not just TEMs)
CD8all<-subset(tcell, subset = final_immune == c("CD8 TEM", "CD8 TCM", "CD8 Proliferating", "CD8 Naive"))
Idents(CD8all)<-"seurat_clusters"
pdf(str_glue("{fig_dir}_{sample_name}_heatmap_exhaustionmarkers_CD8All.pdf"), width=11.6, height=8.2)
features<-c("LAG3", "BTLA", "CTLA4","PDCD1", "HAVCR2")
DoHeatmap(CD8all, features = features, size =3)
dev.off()

Idents(CD8all)<-"seurat_clusters"
pdf(str_glue("{fig_dir}_{sample_name}_heatmap_cytotoxicmarkers_CD8All.pdf"), width=11.6, height=8.2)
features<-c("GNLY", "GZMK", "TNF", "FASLG", "PRF1", "IL2", "IFNG")
DoHeatmap(CD8all, features = features, size =3)
dev.off()

#I think that is the ticket!!
#cluster 6 is mostly CD8 TCM (Ag experienced t cell) expresses a lot of GNLY but very low expression of exhaustion markers
#cluster 3 is mostly CD8 TEM (effector cytotoxic t cells), expresses many cytotoxic markers (mostly GZMK), and low levels of exhaustion markers
#Cluster 1 is mostly CD8 TEM expresses many cytotoxic markers (GNLY, FASLG, PRF1 and INFG), but also expresses a lot of exhaustion markers (LAG3, CTLA4, PDCD1, HAVCR2)

#Also the nebulosa plots are now working.  I will go through Jose's list of plots and see if things are clearer.

#Nebulosa plots
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD4", "FOXP3"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_FOXP3.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()
#Worked
#IFN activated t cells
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "ISG15"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_IFNTCells.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()
#Effector T cells
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "TNFRSF9"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_activatedTCells.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()
#CD127 (IL7R) and tregs
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD4", "IL7R"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_Tregs.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()
#Exhausted T cells
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "CTLA4"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_ExhaustedCD8.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()

p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "HAVCR2"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_ExhaustedCD8HAVCR2.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()

p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "PDCD1"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_ExhaustedCD8PDCD1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()

#More effector markers
p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "GZMB"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_effectorTCellGZMB.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()

p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "PRF1"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_effectorTCellPRF1.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()


p1<- DimPlot(tcell, reduction = "umap", group.by="seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2<-plot_density(tcell, c("CD8A", "GZMK"), joint = TRUE, combine = FALSE)
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_effectorTCellGZMK.pdf"), width=11.6, height=8.2)
plot_grid(p1,p2[[length(p2)]])
dev.off()

#Split_by function to make sure this is all okay.

HN021219A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN021219A")
HN070219A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN070219A")
HN120520A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN120520A")
HN120819A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN120819A")
HN170419A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
HN200519A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN200519A")
HN230620A<-subset(tcell, subset = MajoritySinglet_Individual_Assignment == "HN230620A")

#HN021219A
sample_name<-"HN021219A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN021219A, c("CD8A", "CTLA4"))
dev.off()
#HN070219A
sample_name<-"HN070219A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN070219A, c("CD8A", "CTLA4"))
dev.off()


#HN120520A
sample_name<-"HN120520A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN120520A, c("CD8A", "CTLA4"))
dev.off()
#HN120819A
sample_name<-"HN120819A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN120819A, c("CD8A", "CTLA4"))
dev.off()

#HN170419A
sample_name<-"HN170419A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN170419A, c("CD8A", "CTLA4"))
dev.off()

#HN200519A
sample_name<-"HN200519A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN200519A, c("CD8A", "CTLA4"))
dev.off()

#HN230620A
sample_name<-"HN230620A"
pdf(str_glue("{fig_dir}_{sample_name}_Nebulosa_CTLA4.pdf"), width=11.6, height=8.2)
plot_density(HN230620A, c("CD8A", "CTLA4"))
dev.off()
#Okay, I think I'm ready to rename the t cells.  I'll send an email to Jose so he can look at it, but I'm pretty happy with it.
#cluster 6 is mostly CD8 TCM (Ag experienced t cell) expresses a lot of GNLY but very low expression of exhaustion markers
#cluster 3 is mostly CD8 TEM (effector cytotoxic t cells), expresses many cytotoxic markers (mostly GZMK), and low levels of exhaustion markers
#Cluster 1 is mostly CD8 TEM expresses many cytotoxic markers (GNLY, FASLG, PRF1 and INFG), but also expresses a lot of exhaustion markers (LAG3, CTLA4, PDCD1, HAVCR2)
#I will rename cluster 1 cells which are also CD8 TEM as "CD8 TEM exhausted"

tcell$final_immune<-ifelse((tcell$seurat_clusters == "1" & tcell$final_immune =="CD8 TEM"), "CD8 TEM Exhausted", as.character(tcell$final_immune))

#Worked

#saveout
saveRDS(tcell, file = "integrated_tcell_activated_tcell.rds")

#Create csv files for all the pools
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated"
LN14<-subset(tcell, subset = orig.ident =="LN14")
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

     LN14  LN23rep1  LN23rep2    pri234     pri34 RZ726_LN2 RZ726_PRI   RZ726LN 
     6199      3192      3362      2277      2278       580      7699      6321 

#LN23rep1
LN23rep1<-subset(tcell, subset = orig.ident =="LN23rep1")
sample_name<-"LN23rep1"
annotation<-LN23rep1$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#LN23rep2
LN23rep2<-subset(tcell, subset = orig.ident =="LN23rep2")
sample_name<-"LN23rep2"
annotation<-LN23rep2$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#pri234
pri234<-subset(tcell, subset = orig.ident =="pri234")
sample_name<-"pri234"
annotation<-pri234$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#pri34
pri34<-subset(tcell, subset = orig.ident =="pri34")
sample_name<-"pri34"
annotation<-pri34$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#RZ726_LN2
RZ726_LN2<-subset(tcell, subset = orig.ident =="RZ726_LN2")
sample_name<-"RZ726_LN2"
annotation<-RZ726_LN2$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#RZ726_PRI
RZ726_PRI<-subset(tcell, subset = orig.ident =="RZ726_PRI")
sample_name<-"RZ726_PRI"
annotation<-RZ726_PRI$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#RZ726LN
RZ726LN<-subset(tcell, subset = orig.ident =="RZ726LN")
sample_name<-"RZ726LN"
annotation<-RZ726LN$final_immune
df<-data.frame(annotation)
rownames(df)<-gsub('.{2}$', '', rownames(df))
write.csv(df, file = str_glue("{data_dir}_{sample_name}_tcell_annotation.csv"))

#Now put back into pooled object so can THEN go into master object...I feel like crying today, this is so tedious
#LN14
LN14  LN23rep1  LN23rep2    pri234     pri34 RZ726_LN2 RZ726_PRI   RZ726LN 
     6199      3192      3362      2277      2278       580      7699      6321 
setwd("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_filter_SCT_PCA_simpleannotation")
LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_LN14_tcell_annotation.csv")
rownames(LN14)<-LN14$X
#Read in pooled object
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_filter_SCT_PCA_simpleannotation/LN14_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#save
saveRDS(object, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN14_filter_SCT_PCA_simpleannotation/LN14_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")
#Pull out the .csv files for each individual
sample_name<-"LN14"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN200519A")

#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)

#save
Ind_name<-"HN170419A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN200519A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))


#LN23,rep1
LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_LN23rep1_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep1_filter_SCT_PCA_simpleannotation/LN23_rep1_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep1_filter_SCT_PCA_simpleannotation/LN23_rep1_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")
#Pull out the .csv files for each individual
sample_name<-"LN23_rep1"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN070219A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN120819A")

#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)

#save
Ind_name<-"HN070219A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN120819A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))









#LN23rep2
LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_LN23rep2_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep2_filter_SCT_PCA_simpleannotation/LN23_rep2_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/LN23_rep2_filter_SCT_PCA_simpleannotation/LN23_rep2_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")
#Pull out the .csv files for each individual
sample_name<-"LN23_rep2"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN070219A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN120819A")

#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)

#save
Ind_name<-"HN070219A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN120819A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))


#pri234

LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_pri234_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri234_filter_SCT_PCA_simpleannotation/pri234_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri234_filter_SCT_PCA_simpleannotation/pri234_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")

#Pull out the .csv files for each individual
sample_name<-"pri_234"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN120819A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
p3<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN200519A")
#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)

annotation<-p3$final_immune
df3<-data.frame(annotation)

#save
Ind_name<-"HN120819A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN170419A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN200519A"
write.csv(df3, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))







#pri34

LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_pri34_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri34_filter_SCT_PCA_simpleannotation/pri34_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/pri34_filter_SCT_PCA_simpleannotation/pri34_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")
#Pull out the .csv files for each individual
sample_name<-"pri_34"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN120819A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN200519A")

#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)

#save
Ind_name<-"HN120819A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN200519A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))










#RZ726_LN2

LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_RZ726_LN2_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ726_LN2_filter_SCT_PCA_simpleannotation/RZ726_LN2_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZ726_LN2_filter_SCT_PCA_simpleannotation/RZ726_LN2_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")
#Pull out the .csv files for each individual
sample_name<-"RZ_LN2"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN021219A")


#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)



#save
Ind_name<-"HN021219A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))


#RZ726_PRI

LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_RZ726_PRI_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZpri_filter_SCT_PCA_simpleannotation/RZpri_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZpri_filter_SCT_PCA_simpleannotation/RZpri_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")

#Pull out the .csv files for each individual
sample_name<-"RZ_pri"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN021219A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN120520A")
p3<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN230620A")
#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)


annotation<-p3$final_immune
df3<-data.frame(annotation)

#save
Ind_name<-"HN021219A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN120520A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN230620A"
write.csv(df3, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))



#RZ726LN 
LN14<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/integrated/integrated_RZ726LN_tcell_annotation.csv")
rownames(LN14)<-LN14$X
object<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZLN_filter_SCT_PCA_simpleannotation/RZLN_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation.rds")
object<-AddMetaData(object = object, metadata = LN14)
#Merge
object$final_immune<-ifelse((is.na(object$annotation)), as.character(object$final_immune), as.character(object$annotation))
#worked
object$X<-NULL
object$annotation<-NULL
#saveout
saveRDS(object, file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data/RZLN_filter_SCT_PCA_simpleannotation/RZLN_filter_SCT_PCA_simpleannotation_macrophage_immune_azimuth_annotation_tcell.rds")

#Pull out the .csv files for each individual
sample_name<-"RZ_LN"
table(object$MajoritySinglet_Individual_Assignment)
p1<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN120520A")
p2<-subset(object, subset = MajoritySinglet_Individual_Assignment == "HN230620A")

#Create the .csv files
annotation<-p1$final_immune
df1<-data.frame(annotation)

annotation<-p2$final_immune
df2<-data.frame(annotation)


#save
Ind_name<-"HN120520A"
write.csv(df1, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))
Ind_name<-"HN230620A"
write.csv(df2, file = str_glue("{data_dir}/{sample_name}_{Ind_name}_immunecells_annotation.csv"))








