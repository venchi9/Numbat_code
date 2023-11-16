#Do the cnv bulk clones figures in protrait for the paper
fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate r_4

library(tidyverse)
library(ggplot2)
library(numbat)
library(dplyr)
library(glue)
library(data.table)
library(ggtree)
library(stringr)
library(tidygraph)
library(patchwork)
library(glue)

library(Seurat)
library(viridis)
library(patchwork)
library(clustree)

#These are the samples we are interested in:


sample.name<-"HN120819A"
sample.name<-"HN021219A"
sample.name<-"HN230620A"

#Load in the data
nb <- readRDS(str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}/numbat/numbat.rds"))
#Special one for HN200519A
sample.name<-"HN200519A"
nb<-readRDS("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN200519A/numbat_400/numbat.rds")


#Then print the figure
#Set figure directory
fig_dir<-(str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}/numbat"))
#Special for HN200519A
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN200519A/numbat_400"


pdf(str_glue("{fig_dir}/{sample.name}_cnv_bulk_clones_portrait.pdf", width=8.2, height=11.6))
nb$bulk_clones %>% 
  plot_bulks(
    min_LLR = 50, # filtering CNVs by evidence
    legend = FALSE
  )

  dev.off()

