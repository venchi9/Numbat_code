#Jose will teach me about scPred.
#Currently in this directory:


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

#I spoke with Jose and he said just to use the Azimuth seurat immune cell subtyping.  He said to do it on the pools rather than the individuals 
#This will make subtyping the smaller cell populations easier (like DCs)
#Will need to then manually put in activated T-cells and macrophages
#This is where the simply annotated pools are:

#Directories
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/SimpleAnnotation"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/code"
output<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Azimuth_annotation/data"


#Read in the pool names:
samples <- c("LN14_filter_SCT_PCA_simpleannotation", "LN23_rep1_filter_SCT_PCA_simpleannotation", "LN23_rep2_filter_SCT_PCA_simpleannotation", "pri34_filter_SCT_PCA_simpleannotation", "pri234_filter_SCT_PCA_simpleannotation", "RZ726_LN2_filter_SCT_PCA_simpleannotation", "RZLN_filter_SCT_PCA_simpleannotation","RZpri_filter_SCT_PCA_simpleannotation")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("#Create figure directories")
if(!dir.exists(str_glue('{fig_dir}/{sample_name}'))){
  dir.create(str_glue('{fig_dir}/{sample_name}'))
}

print("#Create output directories")
if(!dir.exists(str_glue('{output}/{sample_name}'))){
  dir.create(str_glue('{output}/{sample_name}'))
}

print("#Load data")
print("load data")
HN<-readRDS(str_glue('{data_dir}/{sample_name}.rds'))

print("subset the immune cells")
print("create column for malignant or non-malignant")
HN$numbat_annotation<-ifelse(((as.character(HN$simple_annotation)=="Cancer") | (as.character(HN$simple_annotation)=="Cancer1")| (as.character(HN$simple_annotation)=="Cancer2")| (as.character(HN$simple_annotation)=="Cancer3")| (as.character(HN$simple_annotation)=="Cancer4")| (as.character(HN$simple_annotation)=="Cancer5")| (as.character(HN$simple_annotation)=="Cancer6")), "Malignant", as.character(HN$simple_annotation))
print("single out stroma and endothelial cells")
HN$numbat_annotation<-ifelse((as.character(HN$numbat_annotation) == "Stroma") | (as.character(HN$numbat_annotation) == "Endothelial") |(as.character(HN$numbat_annotation) == "Malignant"), as.character(HN$numbat_annotation), "Immune")

print("create a subgroup")
immune<-subset(HN, subset = numbat_annotation == "Immune")

print("load reference dataset")
reference<-LoadH5Seurat(str_glue("{data_dir}/pbmc_multimodal.h5seurat"))

print("Find anchors")
anchors <- FindTransferAnchors(
  reference = reference,
  query = immune,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

print("transfer cell type labels")
immune <- MapQuery(
  anchorset = anchors,
  query = immune,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

print("plot the new cell types")
p1 = DimPlot(immune, reduction = "umap", group.by = "simple_annotation", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(immune, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
pdf(str_glue("{fig_dir}/{sample_name}/azimuth_cell_types.pdf"), width=11.6, height=8.2)
p1 + p2
dev.off()

print("save out")
saveRDS(immune, file = str_glue("{output}/{sample_name}/{sample_name}immune_azimuth_annotation.rds"))

print("finished")

