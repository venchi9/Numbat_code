# This is using the new demuliplexed data that Drew and Himanshi did for me in August 2022
#Himanshi put them here /directflow/SCCGGroupShare/projects/himaro/demuxafy/head_neck_venessa/output/combined_results_noSOUP
#I have copied them here:  /directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi
#Note there was mislabeled samples and so the pools are slightly different to the previous experiments
#LN_HN_1_4 - HN170419A and HN200519A **different from before
#LN_HN_2_3 - HN120819A and HN070219A **different from before
#pri_234 - HN120819A, HN170419A, HN200519A (same)
#pri34 - HN120819A, HN200519A (same)
#RZ726_LN - HN120520A, HN230620A (same)
#RZ726_pri - HN021219A, HN120520A, HN230620A
#RZ726_LN2 - HN021219A

#Location of data:
LN_HN_1_4		    /directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_1_4
LN_HN_2_3_rep1		/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_2_3_rep1
LN_HN_2_3_rep2		/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_2_3_rep2
pri_HN_2_3_4		/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/pri_HN_2_3_4
pri_HN_3_4		    /directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/pri_HN_3_4
RZ726_LN		    /directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/200911_VH00259_4_AAAFKY5M5/GE/RZ726_LN
RZ726_PRI	        /directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/200911_VH00259_4_AAAFKY5M5/GE/RZ726_PRI
RZ726_LN2	        /directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/201016_A00152_0317_AHMTN3DSXY/GE/RZ729_LN2

#Have just found out that there was a subdirectory in which the data was put combining the multiple sequencing runs....OMG

/directflow/GWCCGPipeline/projects/bioinformatics/R_191202_VENCHI_INT_10X/analyses/pri_HN_2_3_4/GE/
/directflow/GWCCGPipeline/projects/bioinformatics/R_191202_VENCHI_INT_10X/analyses/LN_HN_1_4/GE/LN_HN_1_4/outs/filtered_feature_bc_matrix
/directflow/GWCCGPipeline/projects/bioinformatics/R_191202_VENCHI_INT_10X/analyses/pri_HN_3_4/GE/
/directflow/GWCCGPipeline/projects/bioinformatics/R_191202_VENCHI_INT_10X/analyses/LN_HN_2_3_rep1/GE/
/directflow/GWCCGPipeline/projects/bioinformatics/R_191202_VENCHI_INT_10X/analyses/LN_HN_2_3_rep2/GE/
/directflow/GWCCGPipeline/projects/bioinformatics/R_200714_VENCHI_INT_10X/analyses/RZ726_PRI/GE/
/directflow/GWCCGPipeline/projects/bioinformatics/R_200714_VENCHI_INT_10X/analyses/RZ726_LN/GE/


#Firstly, need to create Seurat objects and store them.
#Put new objects here:  /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate r_4
cd /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data
R

library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratData)
library(cowplot)
library(Nebulosa)
library(glue)

#Start with LN14
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_1_4/outs/test/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "LN14", min.cells =3, min.features = 200)

#Something wrong with this data:  I have an old raw object here, but maybe this is why the demultiplexing has been hard for this one:
/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/LNintegrated


#LN23_rep1
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_2_3_rep1/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "LN23_rep1", min.cells =3, min.features = 200)
saveRDS(HN_seurat, file = "LN23_rep1_raw.rds")

#LN23_rep2
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_2_3_rep2/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "LN23_rep2", min.cells =3, min.features = 200)
saveRDS(HN_seurat, file = "LN23_rep2_raw.rds")

#pri234
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/pri_HN_2_3_4/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "pri234", min.cells =3, min.features = 200)
saveRDS(HN_seurat, file = "pri234_raw.rds")

#pri34
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/pri_HN_3_4/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "pri34", min.cells =3, min.features = 200)
saveRDS(HN_seurat, file = "pri34_raw.rds")

#RZ726LN
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/200911_VH00259_4_AAAFKY5M5/GE/RZ726_LN/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "RZ726_LN", min.cells =3, min.features = 200)
saveRDS(HN_seurat, file = "RZ726_LN_raw.rds")

#RZ726_PRI **This also doesn't work
HN<- Read10X(data.dir = "/directflow/GWCCGPipeline/projects/bioinformatics/R_200714_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/RZ726_PRI/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "RZ726_PRI", min.cells =3)
saveRDS(HN_seurat, file = "RZ726_PRI_raw.rds")

#RZ726_LN2
HN<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/201016_A00152_0317_AHMTN3DSXY/GE/RZ729_LN2/outs/filtered_feature_bc_matrix")
HN_seurat<-CreateSeuratObject(counts = HN, project = "RZ726_LN2", min.cells =3, min.features = 200)
saveRDS(HN_seurat, file = "RZ726_LN2_raw.rds")

#Now do the QC stuff
#LN23_rep1

pbmc<-readRDS(file = "LN23_rep1_raw.rds")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf("LN23_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#Add in the data from demux.
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/LN_HN_2_3_rep1/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
#Need to set the rownames:
rownames(x)<-x$Barcode

#Add in meta-data
pbmc<-AddMetaData(object = pbmc, metadata = x)
saveRDS(pbmc, file = "LN23_rep1_addfreemux.rds")

************************************************************
#Compare the objects generated from the "analses" folder from Michael

#LN14
old<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_1_4/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "LN14", min.cells =3)
#cells = 17750, max features = 34

new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/analyses/LN_HN_1_4/GE/LN_HN_1_4/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "LN14", min.cells =3)
#cells = 20,000, features = 11004
saveRDS(new_seurat, file = "LN14_raw.rds")

#LN23rep1 
old<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_2_3_rep1/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "LN23rep1", min.cells =3)
#Cells 20,000, max features 7407
new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/analyses/LN_HN_2_3_rep1/GE/LN_HN_2_3_rep1/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "LN23rep1", min.cells =3)
#Cells 20,000, features 9562
saveRDS(new_seurat, file = "LN23_rep1_raw.rds")

#LN23rep2
old<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/LN_HN_2_3_rep2/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "LN23rep2", min.cells =3)
#cells, 20,000, features 7318
new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/analyses/LN_HN_2_3_rep2/GE/LN_HN_2_3_rep2/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "LN23rep2", min.cells =3)
#cells, 20,000, features 10106
saveRDS(new_seurat, file = "LN23_rep2_raw.rds")
#pri234
old<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/pri_HN_2_3_4/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "pri234", min.cells =3)
#Cells - 15,000, max features 7948
new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/analyses/pri_HN_2_3_4/GE/pri_HN_2_3_4/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "pri234", min.cells =3)
#Cells - 15000, max features 10,003
saveRDS(new_seurat, file = "pri234_raw.rds")

#Pri34
old<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/pri_HN_3_4/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "pri34", min.cells =3)
#Cells 20,000, Max features 8899
new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_191202_VENCHI_INT_10X/analyses/pri_HN_3_4/GE/pri_HN_3_4/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "pri34", min.cells =3)
#cells - 20,000 max features 10481
saveRDS(new_seurat, file = "pri34_raw.rds")


#RZ726LN
old<- Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/200911_VH00259_4_AAAFKY5M5/GE/RZ726_LN/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "RZ726LN", min.cells =3)
#Cells 20,000, max features 6849
new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/analyses/RZ726_LN/GE/RZ726_LN/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "RZ726LN", min.cells =3)
#Cells - 20,000, max features 10662
saveRDS(new_seurat,file = "RZ726_LN_raw.rds")

#RZ726Pri
old<-Read10X(data.dir = "/directflow/GWCCGPipeline/projects/bioinformatics/R_200714_VENCHI_INT_10X/200921_VH00259_5_AAAGNV2M5/GE/RZ726_PRI/outs/filtered_feature_bc_matrix")
old_seurat<-CreateSeuratObject(counts = old, project = "RZ726_PRI", min.cells =3)
#Cells 1371, 79 features
new<-Read10X(data.dir = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/R_200714_VENCHI_INT_10X/analyses/RZ726_PRI/GE/RZ726_PRI/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "RZ726_PRI", min.cells =3)
#Cells 20,000, 23538 features
saveRDS(new_seurat, file = "RZ726_PRI_raw.rds")

#Rz726LN2
new<-Read10X(data.dir = "/directflow/GWCCGPipeline/projects/bioinformatics/R_200714_VENCHI_INT_10X/201016_A00152_0317_AHMTN3DSXY/GE/RZ729_LN2/outs/filtered_feature_bc_matrix")
new_seurat<-CreateSeuratObject(counts = new, project = "RZ726_LN2", min.cells =3)
#Cells 20,000, 23538 features
saveRDS(new_seurat, file = "RZ726_LN2_raw.rds")
#Himanshi used the right objects!  Yay...a happy mistake.


#Add in the meta data now:
/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi
#LN14
LN14<-readRDS(file = "LN14_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/LN_HN_1_4/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
LN14<-AddMetaData(object = LN14, metadata = x)
#Add MT 
LN14[["percent.mt"]] <- PercentageFeatureSet(LN14, pattern = "^MT-")
#Plot QC
pdf("LN14_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(LN14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(LN14, file = "LN14_freemux_MT.rds")

#LN23rep1
LN23_1<-readRDS(file = "LN23_rep1_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/LN_HN_2_3_rep1/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
LN23_1<-AddMetaData(object = LN23_1, metadata = x)
#Add MT 
LN23_1[["percent.mt"]] <- PercentageFeatureSet(LN23_1, pattern = "^MT-")
#Plot QC
pdf("LN23_1_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(LN23_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(LN23_1, file = "LN23_rep1_freemux_MT.rds")

#LN23_rep2
LN23_2<-readRDS(file = "LN23_rep2_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/LN_HN_2_3_rep2/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
LN23_2<-AddMetaData(object = LN23_2, metadata = x)
#Add MT 
LN23_2[["percent.mt"]] <- PercentageFeatureSet(LN23_2, pattern = "^MT-")
#Plot QC
pdf("LN23_2_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(LN23_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(LN23_2, file = "LN23_rep2_freemux_MT.rds")


#Pri234
pri234<-readRDS(file = "pri234_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/pri_HN_2_3_4/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
pri234<-AddMetaData(object = pri234, metadata = x)
#Add MT 
pri234[["percent.mt"]] <- PercentageFeatureSet(pri234, pattern = "^MT-")
#Plot QC
pdf("pri234_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(pri234, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(pri234, file = "pri234_freemux_MT.rds")

#pri34
pri34<-readRDS(file = "pri34_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/pri_HN_3_4/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
pri34<-AddMetaData(object = pri34, metadata = x)
#Add MT 
pri34[["percent.mt"]] <- PercentageFeatureSet(pri34, pattern = "^MT-")
#Plot QC
pdf("pri34_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(pri34, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(pri34, file = "pri34_freemux_MT.rds")

#RZ726LN
RZLN<-readRDS(file = "RZ726_LN_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/RZ726_LN/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
RZLN<-AddMetaData(object = RZLN, metadata = x)
#Add MT 
RZLN[["percent.mt"]] <- PercentageFeatureSet(RZLN, pattern = "^MT-")
#Plot QC
pdf("RZ726_LN_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(RZLN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(RZLN, file = "RZ726_LN_freemux_MT.rds")

#RZ726_PRI
RZpri<-readRDS(file = "RZ726_PRI_raw.rds")
x<-read.table(file = "/directflow/SCCGGroupShare/projects/venchi/HN/Raw_Data/HN_demultiplexing_July2022_Himanshi/RZ726_PRI/combined_results_w_combined_assignments.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
RZpri<-AddMetaData(object = RZpri, metadata = x)
#Add MT 
RZpri[["percent.mt"]] <- PercentageFeatureSet(RZpri, pattern = "^MT-")
#Plot QC
pdf("RZ726_pri_rep1_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(RZpri, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(RZpri, file = "RZ726_pri_freemux_MT.rds")

#RZ726LN2
x<-read.table(file = "/directflow/SCCGGroupShare/projects/himaro/demuxafy/head_neck_venessa/output/scds/RZ729_LN2/scds_doublets_singlets.tsv", sep = '\t', header = TRUE)
rownames(x)<-x$Barcode
new_seurat<-AddMetaData(object = new_seurat, metadata = x)
#Add MT 
new_seurat[["percent.mt"]] <- PercentageFeatureSet(new_seurat, pattern = "^MT-")
#Plot QC
pdf("RZ726_LN2_QC_no_filtering.pdf", width=11.6, height=8.2)
VlnPlot(new_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Save
saveRDS(new_seurat, file = "RZ726_LN2_freemux_MT.rds")

#Now I have to add viral reads and filter.
#What a pain.

# I have previously extracted the reads into these text files:  
/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/

#LN14
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/LN14"

# Import HPV-mapped UMI info
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
LN14$hpv_all <- hpv$n[match(str_split_fixed(colnames(LN14), "_", 2)[,1], hpv$CB)]
LN14$hpv_all[is.na(LN14$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
LN14<-subset(LN14, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
LN14<-subset(LN14, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
LN14<-subset(LN14, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#QC plot
pdf("LN14_rep1_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(LN14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(LN14, file = "LN14_freemux_MT_filtering_viralreads.rds")

#LN23rep1
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/LN23rep1"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
LN23_1$hpv_all <- hpv$n[match(str_split_fixed(colnames(LN23_1), "_", 2)[,1], hpv$CB)]
LN23_1$hpv_all[is.na(LN23_1$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
LN23_1<-subset(LN23_1, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
LN23_1<-subset(LN23_1, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
LN23_1<-subset(LN23_1, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#Change "donor 1" to HN070219A
LN23_1$MajoritySinglet_Individual_Assignment<-ifelse(as.character(LN23_1$MajoritySinglet_Individual_Assignment == "donor1"), "HN070219A", as.character(LN23_1$MajoritySinglet_Individual_Assignment))
#Save
pdf("LN23_rep1_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(LN23_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(LN23_1, file = "LN23_rep1_freemux_MT_filtering_viralreads.rds")

#LN23_rep2
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/LN23rep2"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
LN23_2$hpv_all <- hpv$n[match(str_split_fixed(colnames(LN23_2), "_", 2)[,1], hpv$CB)]
LN23_2$hpv_all[is.na(LN23_2$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
LN23_2<-subset(LN23_2, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
LN23_2<-subset(LN23_2, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
LN23_2<-subset(LN23_2, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#Change "donor 1" to HN070219A
LN23_2$MajoritySinglet_Individual_Assignment<-ifelse(as.character(LN23_2$MajoritySinglet_Individual_Assignment == "donor1"), "HN070219A", as.character(LN23_2$MajoritySinglet_Individual_Assignment))
#Save
pdf("LN23_2_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(LN23_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(LN23_2, file = "LN23_rep2_freemux_MT_filtering_viralreads.rds")

#PRi234
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/pri234"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
pri234$hpv_all <- hpv$n[match(str_split_fixed(colnames(pri234), "_", 2)[,1], hpv$CB)]
pri234$hpv_all[is.na(pri234$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
pri234<-subset(pri234, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
pri234<-subset(pri234, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
pri234<-subset(pri234, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#Save
pdf("pri234_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(pri234, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(pri234, file = "pri234_freemux_MT_filtering_viralreads.rds")

#pri34
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/pri34"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
pri34$hpv_all <- hpv$n[match(str_split_fixed(colnames(pri34), "_", 2)[,1], hpv$CB)]
pri34$hpv_all[is.na(pri34$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
pri34<-subset(pri34, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
pri34<-subset(pri34, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
pri34<-subset(pri34, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#Save
pdf("pri34_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(pri34, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(pri34, file = "pri34_freemux_MT_filtering_viralreads.rds")

#RZ726_LN
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/RZ726LN"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
RZLN$hpv_all <- hpv$n[match(str_split_fixed(colnames(RZLN), "_", 2)[,1], hpv$CB)]
RZLN$hpv_all[is.na(RZLN$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
RZLN<-subset(RZLN, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
RZLN<-subset(RZLN, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
RZLN<-subset(RZLN, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#Save
pdf("RZLN_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(RZLN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(RZLN, file = "RZLN_freemux_MT_filtering_viralreads.rds")

#RZ726_pri
output_dir <- "//directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/RZ726P"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
RZpri$hpv_all <- hpv$n[match(str_split_fixed(colnames(RZpri), "_", 2)[,1], hpv$CB)]
RZpri$hpv_all[is.na(RZpri$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
RZpri<-subset(RZpri, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
RZpri<-subset(RZpri, subset = MajoritySinglet_Individual_Assignment != "doublet")
#Remove unassigned
RZpri<-subset(RZpri, subset = MajoritySinglet_Individual_Assignment != "unassigned")
#Save
pdf("RZpri_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(RZpri, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(RZpri, file = "RZpri_freemux_MT_filtering_viralreads.rds")

#RZ726_LN2
output_dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HNSinglecell/FinalAnalysisOct2020/Addviralreads/RZ726LN2"
hpv <- read_tsv(glue('{output_dir}/CB_UB.txt'), col_names = c("CB", "UB")) %>%
  group_by(CB) %>%          
  summarise(n = n_distinct(UB)) 
head(hpv)

# Add to Seurat object
new_seurat$hpv_all <- hpv$n[match(str_split_fixed(colnames(new_seurat), "_", 2)[,1], hpv$CB)]
new_seurat$hpv_all[is.na(new_seurat$hpv_all)] <- 0

#Now do filtering with these paramters:  nFeatures>200, nCounts >200, MT<25%
new_seurat<-subset(new_seurat, subset = nFeature_RNA >200 & nCount_RNA>200 & percent.mt<25)

#Remove doublets
new_seurat<-subset(new_seurat, subset = scds_DropletType != "doublet")

#Create MajoritySinglet_Individual_Assignment
new_seurat$MajoritySinglet_Individual_Assignment<-"HN021219A"

#Save
pdf("RZ726_LN2_QC_filtering.pdf", width=11.6, height=8.2)
VlnPlot(new_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "hpv_all"), ncol = 4)
dev.off()
saveRDS(new_seurat, file = "RZ726_LN2_freemux_MT_filtering_viralreads.rds")





#Generate list of barcodes for each patient for each object
#LN14
LN14_2005<-subset(LN14, subset = MajoritySinglet_Individual_Assignment =="HN200519A")
names<-colnames(LN14_2005)
write.csv(names, file = "LN14_HN200519A_barcodes.csv")
LN14_1704<-subset(LN14, subset = MajoritySinglet_Individual_Assignment =="HN170419A")
names<-colnames(LN14_1704)
write.csv(names, file = "LN14_HN170419A_barcodes.csv")

#LN23rep1
LN231_1208<-subset(LN23_1, subset = MajoritySinglet_Individual_Assignment =="HN120819A")
names<-colnames(LN231_1208)
write.csv(names, file = "LN23_rep1_HN120819A")
LN231_0702<-subset(LN23_1, subset = MajoritySinglet_Individual_Assignment =="HN070219A")
names<-colnames(LN231_0702)
write.csv(names, file = "LN23_rep1_HN070219A.csv")

#LN23rep2
LN232_1208<-subset(LN23_2, subset = MajoritySinglet_Individual_Assignment =="HN120819A")
names<-colnames(LN232_1208)
write.csv(names, file = "LN23_rep2_HN120819A")
LN232_0702<-subset(LN23_2, subset = MajoritySinglet_Individual_Assignment =="HN070219A")
names<-colnames(LN232_0702)
write.csv(names, file = "LN23_rep2_HN070219A.csv")

#pri234
pri234_1208<-subset(pri234, subset = MajoritySinglet_Individual_Assignment =="HN120819A")
names<-colnames(pri234_1208)
write.csv(names, file = "pri234_HN120819A.csv")

pri234_1208<-subset(pri234, subset = MajoritySinglet_Individual_Assignment =="HN170419A")
names<-colnames(pri234_1208)
write.csv(names, file = "pri234_HN170419A.csv")

pri234_1208<-subset(pri234, subset = MajoritySinglet_Individual_Assignment =="HN200519A")
names<-colnames(pri234_1208)
write.csv(names, file = "pri234_HN200519A.csv")

#pri34
pri34_1208<-subset(pri34, subset = MajoritySinglet_Individual_Assignment =="HN200519A")
names<-colnames(pri34_1208)
write.csv(names, file = "pri34_HN200519A.csv")

pri34_1208<-subset(pri34, subset = MajoritySinglet_Individual_Assignment =="HN120819A")
names<-colnames(pri34_1208)
write.csv(names, file = "pri34_HN120819A.csv")

#RZLN
RZLN_23<-subset(RZLN, subset = MajoritySinglet_Individual_Assignment =="HN230620A")
names<-colnames(RZLN_23)
write.csv(names, file = "RZ726_LN_HN230620A.csv")

RZLN_23<-subset(RZLN, subset = MajoritySinglet_Individual_Assignment =="HN120520A")
names<-colnames(RZLN_23)
write.csv(names, file = "RZ726_LN_HN120520A.csv")

#RZpri
RZpri_23<-subset(RZpri, subset = MajoritySinglet_Individual_Assignment =="HN021219A")
names<-colnames(RZpri_23)
write.csv(names, file = "RZ726_pri_HN021219A.csv")

RZpri_23<-subset(RZpri, subset = MajoritySinglet_Individual_Assignment =="HN120520A")
names<-colnames(RZpri_23)
write.csv(names, file = "RZ726_pri_HN120520A.csv")

RZpri_23<-subset(RZpri, subset = MajoritySinglet_Individual_Assignment =="HN230620A")
names<-colnames(RZpri_23)
write.csv(names, file = "RZ726_pri_HN230620A.csv")

#RZ726_LN2
names<-colnames(new_seurat)
write.csv(names, file = "RZ726_LN2_HN021219A.csv")

new_seurat <- SCTransform(new_seurat, method = "glmGamPoi", verbose = FALSE)
new_seurat
head(new_seurat@meta.data)

#Now I can annotate the cells roughly using the feature plots that I generated using a shell code
#Objects are still here:/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data

#start with LN14
LN14<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/LN14_filter_SCT_PCA.rds")
new.cluster.ids<-c("NK", "CD8_1", "CD4_1", "CD4_2", "BCells", "Cancer1", "CD8_2", "Cancer2", "Cancer3", "Cancer4", "Mono", "CD8_3")
names(new.cluster.ids) <- levels(LN14)
LN14 <- RenameIdents(LN14, new.cluster.ids)
LN14[["simple_annotation"]] <- Idents(LN14)
#Plot UMAP
pdf("LN14_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(LN14, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend()
dev.off()
#saveout
saveRDS(LN14, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/LN14_filter_SCT_PCA_simpleannotation.rds")

#LN23_rep1
LN23_1<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/LN23_rep1_filter_SCT_PCA.rds")
new.cluster.ids<-c("Cancer1", "Cancer2", "BCells", "CD4_1", "CD8_1", "CD8_2", "CD4_2", "Cancer3", "Cancer4", "CD8_3", "Mono", "Endothelial")
names(new.cluster.ids) <- levels(LN23_1)
LN23_1 <- RenameIdents(LN23_1, new.cluster.ids)
LN23_1[["simple_annotation"]] <- Idents(LN23_1)
#Plot UMAP
pdf("LN23_1_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(LN23_1, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend()
dev.off()

#Saveout
saveRDS(LN23_1, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/LN23_rep1_filter_SCT_PCA_simpleannotation.rds")

#LN23_rep2
LN23_2<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/LN23_rep2_filter_SCT_PCA.rds")
new.cluster.ids<-c("Cancer1", "BCells", "Cancer2", "CD8_1", "CD4_1", "CD4_2", "Cancer3", "CD8_2", "CD4_3", "Cancer4", "Mono", "Endothelial", "Stroma")
names(new.cluster.ids) <- levels(LN23_2)
LN23_2 <- RenameIdents(LN23_2, new.cluster.ids)
LN23_2[["simple_annotation"]] <- Idents(LN23_2)
#Plot UMAP
pdf("LN23_2_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(LN23_2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(LN23_2, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/LN23_rep2_filter_SCT_PCA_simpleannotation.rds")

#pri34
pri34<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/pri34_filter_SCT_PCA.rds")
new.cluster.ids<-c("BCells", "Cancer1", "CD4_1", "Cancer2", "CD4_2", "CD8_1", "Cancer3", "CD8_2", "CD8_3", "CD8_4", "Cancer4", "Mono", "Cancer5", "Cancer6", "Endothelial", "Plasma")
names(new.cluster.ids) <- levels(pri34)
pri34 <- RenameIdents(pri34, new.cluster.ids)
pri34[["simple_annotation"]] <- Idents(pri34)
#Plot UMAP
pdf("pri34_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(pri34, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(pri34, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/pri34_filter_SCT_PCA_simpleannotation.rds")

#pri234
pri234<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/pri234_filter_SCT_PCA.rds")
new.cluster.ids<-c("BCells", "Cancer1", "CD8_1", "CD4_1", "CD4_2", "Cancer2", "CD8_2", "Cancer3", "Cancer4", "Endothelial", "Mono", "Cancer5", "Plasma", "Cancer6", "Stroma")
names(new.cluster.ids) <- levels(pri234)
pri234 <- RenameIdents(pri234, new.cluster.ids)
pri234[["simple_annotation"]] <- Idents(pri234)
#Plot UMAP
pdf("pri234_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(pri234, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(pri234, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/pri234_filter_SCT_PCA_simpleannotation.rds")

#RZ726_LN2

LN2<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/RZ726_LN2_filter_SCT_PCA.rds")
new.cluster.ids<-c("CD4", "Cancer1", "Cancer2", "Cancer3", "Cancer4", "Cancer5", "Stroma", "Endothelial", "CD8", "Cancer6", "Mono", "BCells")
names(new.cluster.ids) <- levels(LN2)
LN2 <- RenameIdents(LN2, new.cluster.ids)
LN2[["simple_annotation"]] <- Idents(LN2)
#Plot UMAP
pdf("RZ726_LN2_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(LN2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(LN2, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/RZ726_LN2_filter_SCT_PCA_simpleannotation.rds")

#RZ726_LN
LN<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/RZLN_filter_SCT_PCA.rds")
new.cluster.ids<-c("BCell1", "CD4_1", "CD8_1", "Cancer1" ,"CD4_2", "CD4_3", "BCell2", "BCell3", "CD8_2", "CD8_3", "Mono", "NK", "DC", "Stroma")
names(new.cluster.ids) <- levels(LN)
LN <- RenameIdents(LN, new.cluster.ids)
LN[["simple_annotation"]] <- Idents(LN)
#Plot UMAP
pdf("RZ726_LN_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(LN, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(LN, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/RZLN_filter_SCT_PCA_simpleannotation.rds")

#RZ726_pri
RZ<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/RZpri_filter_SCT_PCA.rds")
new.cluster.ids<-c("BCell1", "CD4_1", "CD8_1", "CD4_2", "BCell2", "CD8_2", "CD4_3", "Cancer", "Mono", "BCell3", "CD8_3", "Endothelial", "Stroma", "DC", "Unknown")
names(new.cluster.ids) <- levels(RZ)
RZ <- RenameIdents(RZ, new.cluster.ids)
RZ[["simple_annotation"]] <- Idents(RZ)
#Plot UMAP
pdf("RZ726_pri_UMAP_annotated.pdf", width=11.6, height=8.2)
DimPlot(RZ, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
saveRDS(RZ, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/RZpri_filter_SCT_PCA_simpleannotation.rds")

#Now need to split apart the individual patients
#Start with LN14
LN14$Location<-"Lymph_node"
HN1704<-subset(LN14, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
HN200<-subset(LN14, subset = MajoritySinglet_Individual_Assignment == "HN200519A")
#Saveout
saveRDS(HN1704, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN170419A/LN14_HN170419A.rds")
saveRDS(HN200,file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN200519A/LN14_HN200519A.rds")

#LN23rep1
LN23_1$Location<-"Lymph_node"
HN1208<-subset(LN23_1, subset = MajoritySinglet_Individual_Assignment =="HN120819A")
HN0702<-subset(LN23_1, subset = MajoritySinglet_Individual_Assignment =="HN070219A")
#Save out
saveRDS(HN1208, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN120819A/LN23_1_HN120819A.rds")
saveRDS(HN0702, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN070219A/LN23_1_HN070219A.rds")

#LN23_rep2
LN23_2$Location<-"Lymph_node"
HN1208<-subset(LN23_2, subset = MajoritySinglet_Individual_Assignment =="HN120819A")
HN0702<-subset(LN23_2, subset = MajoritySinglet_Individual_Assignment =="HN070219A")
#save
saveRDS(HN1208, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN120819A/LN23_2_HN120819A.rds")
saveRDS(HN0702, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN070219A/LN23_2_HN070219A.rds")

#RZ726_LN
LN$Location<-"Lymph_node"
HN2306<-subset(LN, subset = MajoritySinglet_Individual_Assignment == "HN230620A")
HN1205<-subset(LN, subset = MajoritySinglet_Individual_Assignment == "HN120520A")
#Save
saveRDS(HN2306, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN230620A/RZ726_LN_HN230620A.rds")
saveRDS(HN1205, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN120520A/RZ726_LN_HN120520A.rds")

#RZ726_LN2
LN2$Location<-"Lymph_node"
saveRDS(LN2, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN021219A/RZ726_LN2_HN021219A")

#Pri234
pri234$Location<-"Primary"
HN1208<-subset(pri234, subset = MajoritySinglet_Individual_Assignment == "HN120819A")
HN1704<-subset(pri234, subset = MajoritySinglet_Individual_Assignment == "HN170419A")
HN2005<-subset(pri234, subset = MajoritySinglet_Individual_Assignment == "HN200519A")

#Save
saveRDS(HN1208, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN120819A/pri234_HN120819A")
saveRDS(HN1704, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN170419A/pri234_HN170419A")
saveRDS(HN2005, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN200519A/pri234_HN200519A")

#pri34
pri34$Location<-"Primary"
HN1208<-subset(pri34, subset = MajoritySinglet_Individual_Assignment == "HN120819A")
HN2005<-subset(pri34, subset = MajoritySinglet_Individual_Assignment == "HN200519A")

#Save
saveRDS(HN1208, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN120819A/pri34_HN120819A.rds")
saveRDS(HN2005, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN200519A/pri34_HN200519A.rds")

#RZ726_pri
RZ$Location<-"Primary"
HN0212<-subset(RZ, subset = MajoritySinglet_Individual_Assignment =="HN021219A")
HN1205<-subset(RZ, subset = MajoritySinglet_Individual_Assignment =="HN120520A")
HN2306<-subset(RZ, subset = MajoritySinglet_Individual_Assignment =="HN230620A")

#Save
saveRDS(HN0212,file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN021219A/RZ726_pri_HN021219A.rds")
saveRDS(HN1205, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN120520A/RZ726_pri_HN120520A.rds")
saveRDS(HN2306, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/HN230620A/RZ726_pri_HN230620A.rds")

#I tidied up the data folder a bit.  All the pre-processing objects are now in a folder called "Data"

#Merge all the objects together into individuals.

#Data is here: /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data
#HN021219A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN021219A/RZ726_pri_HN021219A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN021219A/RZ726_LN2_HN021219A.rds")

comb<- list(a, b)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="Location")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN021219A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Save
saveRDS(exp.integrated, file = "HN021219A_integrated_object.rds")

#HN070219A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN070219A/LN23_2_HN070219A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN070219A/LN23_1_HN070219A.rds")

comb<- list(a, b)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="orig.ident")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN070219A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(exp.integrated, file = "HN070219A_integrated_object.rds")

#HN120520A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120520A/RZ726_pri_HN120520A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120520A/RZ726_LN_HN120520A.rds")

comb<- list(a, b)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="Location")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN120520A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(exp.integrated, file = "HN120520A_integrated_object.rds")

#HN120819A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/pri34_HN120819A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/pri234_HN120819A.rds")
c<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/LN23_2_HN120819A.rds")
d<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN120819A/LN23_1_HN120819A.rds")


comb<- list(a, b, c, d)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b, c, d), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b, c, d), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="Location")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN120819A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#Save
saveRDS(exp.integrated, file = "HN120819A_integrated_object.rds")

#HN170419A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN170419A/pri234_HN170419A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN170419A/LN14_HN170419A.rds")

comb<- list(a, b)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="Location")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN170419A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(exp.integrated, file = "HN170419A_integrated_object.rds")

#HN200519A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A/pri34_HN200519A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A/pri234_HN200519A.rds")
c<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN200519A/LN14_HN200519A.rds")

comb<- list(a, b,c)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b,c), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b,c), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="Location")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN200519A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(exp.integrated, file = "HN200519A_integrated_object.rds")

#HN230620A
a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN230620A/RZ726_pri_HN230620A.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/HN230620A/RZ726_LN_HN230620A.rds")

comb<- list(a, b)
Int.feat<- SelectIntegrationFeatures(object.list = list(a, b), nfeatures = 3000)
comb<- PrepSCTIntegration(object.list = list(a, b), anchor.features = Int.feat)
int.anchors <- FindIntegrationAnchors(object.list = comb,normalization.method = "SCT", anchor.features = Int.feat)
exp.integrated<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

#Do PCA and umap
exp.integrated <- RunPCA(object = exp.integrated, npcs = 100, verbose = FALSE)
exp.integrated<-RunUMAP(exp.integrated, reduction = "pca", dims = 1:15)
exp.integrated<-FindNeighbors(exp.integrated, reduction= "pca", dims = 1:15)
exp.integrated <- FindClusters(exp.integrated, resolution = 0.5)

#Plot
Idents(exp.integrated)<-"simple_annotation"
p1<- DimPlot(exp.integrated, reduction = "umap", group.by="Location")
p2<-DimPlot(exp.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
pdf("HN230620A_integrated_UMAP.pdf", width=11.6, height=8.2)
plot_grid(p1,p2)
dev.off()

#save
saveRDS(exp.integrated, file = "HN230620A_integrated_object.rds")

