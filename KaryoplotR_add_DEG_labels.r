#Do KaryoploteR on all the samples for the paper and add the labels of the DEG

suppressPackageStartupMessages({
library(karyoploteR)
library(numbat)
library(tidyverse)
library(biomaRt)
})

# Define directory
data.dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

print("set sample names")
samples <- c("HN021219A", "HN120819A", "HN200519A", "HN230620A", "HN070219A", "combined_cnv")
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])
sample_name <- samples[i]

print("load in the data")
cnv<-read_csv(file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/{sample_name}_CNV_summary.csv"))

print("start plotting")
pp <- getDefaultPlotParams(plot.type = 3)
pp$leftmargin <- 0.15

# Define some genes of interest and get their genomic coordinates using biomaRt
gene.symbols <- c("RPL5", "RPS24", "RPL26")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes_add_DEG_labels.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr{1:22}'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.25, r1 = 0.5)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.8,0.1,0.1),
              r1=0.7)

# Close the plot
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes_mViralRNA_labels.pdf"), width=11.6, height=8.2)

gene.symbols <- c("RPL13" , "RPS2",   "RPL10",  "RPS27A", "RPL30" , "RPL12",  "RPS6" ,  "RPL7A" ,
 "RPL5",   "RPSA" ,  "RPS7",  "RPLP1",  "RPL15",  "RPS24" , "RPL23A", "RPL19", "RPL29" , "RPL26" , "RPL41",  "RPL35A")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr{1:22}'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.25, r1 = 0.5)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.8,0.1,0.1),
              r1=0.7)

# Close the plot
dev.off()

#Do plot with other important markers
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes_canonicalmarkers_labels.pdf"), width=11.6, height=8.2)

gene.symbols <- c("CCND1", "FGF4", "CDKN2A", "FGF19", "CDKN2B", "PIK3CA", "MYC")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr{1:22}'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.25, r1 = 0.5)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.8,0.1,0.1),
              r1=0.7)

# Close the plot
dev.off()

#Do plot with other important markers taken from: https://www.nature.com/articles/nature14129
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes_TCGA_markers.pdf"), width=11.6, height=8.2)

gene.symbols <- c("PIK3CA", "TRAF3", "E2F1")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr{1:22}'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 2,
              col = "#cb181d",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="amplification", data.panel =2, r1 = 0)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0.26, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.25)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.51, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 0.5)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.8,0.3,0.3),
              r1=0.3)

# Close the plot
dev.off()

#Do plot with other important markers taken from: https://www.nature.com/articles/nature14129; but also summarised here:  https://link.springer.com/chapter/10.1007/978-3-031-23175-9_6 
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes_olearymarkers.pdf"), width=11.6, height=8.2)

gene.symbols <- c("PIK3CA", "TRAF3", "E2F1", "MYC", "CDKN2A", "RB1", "TP53", "E542K", "E545K")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr{1:22}'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 2,
              col = "#cb181d",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="amplification", data.panel =2, r1 = 0)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0.26, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.25)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.51, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 0.5)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.5,0.3,0.3),
              r1=0.3)

# Close the plot
dev.off()

#Do plot with other important markers taken from: https://www.nature.com/articles/nature14129; but also summarised here:  https://link.springer.com/chapter/10.1007/978-3-031-23175-9_6, also add in tranlsational markers too:  https://pubmed.ncbi.nlm.nih.gov/34358439/

pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes_olearymarkers.pdf"), width=11.6, height=8.2)

gene.symbols <- c("PIK3CA", "TRAF3", "E2F1", "MYC", "CDKN2A", "RB1", "TP53", "E542K", "E545K", "EIF4G1", "EIF4EBP1", "EIF4A1", "EIF4G2", "EIF4E3", "MKNK1")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr{1:22}'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 2,
              col = "#cb181d",
              r0 = 0, r1 = 0.25)
kpAddLabels(kp, labels="amplification", data.panel =2, r1 = 0)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0.26, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.25)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.51, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 0.5)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.2,0.3,0.3),
              r1=0.3)

# Close the plot
dev.off()