# Script information ------------------------------------------------------

# title: Run Numbat R prep
# author: Walter Muskovic 
# date: 2022-07-26
# description: This script will prepare the files needed to run numbat



# Load R packages ---------------------------------------------------------

library(numbat)
library(Seurat)
library(tidyverse)



# Define samples and directories ------------------------------------------

i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
samples <- c("HN200519A") #HN200519A needs a third run with more stringent parameters #HN021219A done during our test and completed. #"HN070219A","HN120520A","HN120819A","HN170419A","HN230620A" 
#For 200519A and 070219A the clones aren't as stable, so will take those to 35 iterations.
sample.name <- samples[i]
seurat.dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
data.dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data"
print(sample.name)



# Prepare data ------------------------------------------------------------

# Get counts for all cells
hn <- readRDS(str_glue('{seurat.dir}/{sample.name}/{sample.name}_numbat_object.rds'))
# Remove duplicated cell IDs
hn <- subset(hn, cell=Cells(hn)[!duplicated(str_split_fixed(Cells(hn),"_",2)[,1])])
# Remove appended "-2"
hn <- RenameCells(hn, new.names=str_split_fixed(Cells(hn),"_",2)[,1])
print("Cell types present in current sample:")
print(table(hn$numbat_annotation))
hn.counts <- as.matrix(hn@assays$RNA@counts)

# Get vectors with the cell barcodes of all normal and malignant cells
normal.cells <- colnames(hn)[hn$numbat_annotation!="Malignant"]
malignant.cells <- colnames(hn)[hn$numbat_annotation=="Malignant"]

## Make reference gene expression profiles 
print("Creating reference gene expression profiles")
cell.types <- data.frame(cell = colnames(hn),
                         numbat_annotation = hn$numbat_annotation) %>%
    filter(numbat_annotation != "Malignant")

# remove cell types with less than 30 cells
cell.type.counts <- data.frame(table(cell.types$numbat_annotation))
colnames(cell.type.counts) <- c("numbat_annotation", "count")
cell.types <- merge(x = cell.types, y = cell.type.counts, by = "numbat_annotation")
cell.types <- filter(cell.types, count >=30) %>%
    select(-count)
print("Using the following cell types as references:")
print(table(cell.types$numbat_annotation))
# Make reference
colnames(cell.types) <- c("group", "cell")
ref.df <- aggregate_counts(count_mat = hn.counts[,colnames(hn.counts)%in%cell.types$cell],
                           annot = cell.types,
                           normalized = TRUE,
                           verbose = TRUE)

# Get allele counts
print("Importing allele counts")
allele.counts <- read_tsv(str_glue('{data.dir}/{sample.name}/{sample.name}_allele_counts.tsv.gz'))

# Remove normal cells from both the count_matrix and df_allele
print("Removing normal cells from count_matrix and df_allele")
print("Dimensions of counts_matrix with all cells:")
print(dim(hn.counts))
hn.counts <- hn.counts[,colnames(hn.counts)%in%malignant.cells]
print("Dimensions of counts_matrix with just malignant cells:")
print(dim(hn.counts))
print("Dimensions of df_allele with all cells:")
print(dim(allele.counts))
print("Dimensions of df_allele with just malignant cells:")
allele.counts <- allele.counts[allele.counts$cell%in%malignant.cells,]
print(dim(allele.counts))



# Detecting clonal LOH ----------------------------------------------------
bulk <- get_bulk(
    count_mat = hn.counts,
    lambdas_ref = ref.df,
    df_allele = allele.counts,
    gtf = gtf_hg38,
    genetic_map = genetic_map_hg38
)

segs_loh <- bulk %>% detect_clonal_loh(t = 1e-3)

message("Finished preparing data to run numbat")



# Run numbat --------------------------------------------------------------

options(future.globals.maxSize=1e20)

message("Started running numbat")

  out = run_numbat(
   count_mat = hn.counts, # gene x cell integer UMI count matrix 
   lambdas_ref = ref.df, # reference expression profile, a gene x cell type normalized expression level matrix
   df_allele = allele.counts, # allele dataframe generated by pileup_and_phase script
   gtf_hg38, # provided upon loading the package
   genetic_map_hg38, # provided upon loading the package
   min_cells = 20,
   t = 1e-3,
    max_iter = 20,
    min_LLR = 450,
    init_k = 3,
    ncores = 16,
    ncores_nni = 16,
    plot = FALSE,
    out_dir = str_glue('{data.dir}/{sample.name}/numbat_450')
)

message("Finished running numbat")
#Second run - with HN070219A and HN200519A the clones were unstable.  I have changed the min_LLR to 60 and the tau to 0.9
#Original settings here:annotout = run_numbat(
  #  count_mat = hn.counts, # gene x cell integer UMI count matrix 
  #  lambdas_ref = ref.df, # reference expression profile, a gene x cell type normalized expression level matrix
  #  df_allele = allele.counts, # allele dataframe generated by pileup_and_phase script
  #  gtf_hg38, # provided upon loading the package
  #  genetic_map_hg38, # provided upon loading the package
  #  min_cells = 20,
  #  t = 1e-3,
  #  max_iter = 20,
  #  min_LLR = 50,
  #  init_k = 3,
  #  ncores = 16,
  #  ncores_nni = 16,
  #  plot = FALSE,
  #  out_dir = str_glue('{data.dir}/{sample.name}/numbat')
#)

#Third run - HN200519A is still very unstable with the second run parameters below.  I will make them more stringent
#out = run_numbat(
#    count_mat = hn.counts, # gene x cell integer UMI count matrix 
#    lambdas_ref = ref.df, # reference expression profile, a gene x cell type normalized expression level matrix
#    df_allele = allele.counts, # allele dataframe generated by pileup_and_phase script
#    gtf_hg38, # provided upon loading the package
#    genetic_map_hg38, # provided upon loading the package
#    min_cells = 20,
#    t = 1e-3,
#    max_iter = 20,
#    min_LLR = 60,
#    init_k = 3,
#    ncores = 16,
#    ncores_nni = 16,
#    tau = 0.9,
#    plot = FALSE,
#    out_dir = str_glue('{data.dir}/{sample.name}/numbat')
#)

