#try karyoploteR for chromosomal mapping of CNV
# https://bioconductor.org/packages/release/bioc/html/karyoploteR.html
#Walter did this code

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2

library(karyoploteR)
library(numbat)
library(tidyverse)
library(biomaRt)

# Define directory
data.dir <- "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

# Define samples
samples <- c("HN021219A", "HN070219A", "HN120520A", "HN120819A", "HN170419A", "HN230620A", "HN200519A")
#Using numbat_400 for HN200519A

# Load the CNV coordinates from each numbat run
cnv <- lapply(samples, function(sample.name){
  message(sample.name)
  cnv <- readRDS(str_glue("{data.dir}/{sample.name}/numbat/numbat.rds"))$joint_post %>%
    select(CHROM, cnv_state, seg_start, seg_end, seg, LLR) %>%
    distinct() %>%
    mutate(CHROM = str_glue('chr{CHROM}'),
           sample = sample.name)
  return(cnv)
  }) %>%
  do.call(rbind,.)

  #This won't work for me.  Try manually
  sample.name<-"HN021219A"
  sample.name<-"HN200519A"
  sample.name<-"HN120819A"
  sample.name<-"HN230620A"
setwd(str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}/numbat"))
#Load up numbat outs
nb = Numbat$new(str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}/numbat"))
#Create the data.frame
cnv<-nb$joint_post %>% dplyr::select(CHROM, cnv_state, seg_start, seg_end, seg, LLR)%>%
distinct()%>%
mutate(CHROM = str_glue("chr{CHROM}"))
cnv$individual<-sample.name
#save out
write.csv(cnv, file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/{sample.name}_CNV_summary.csv"))

#Now read in and bind together.
my_list<-list.files("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data", full.names = TRUE)
my_files<-lapply(my_list, read_csv)
df<-bind_rows(my_files)
df

write.csv(df, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/combined_cnv.csv")

#Can try walter's fancy way and see if it's the same:
samples <- c("HN021219A",   "HN120819A",  "HN230620A", "HN200519A")

# Load the CNV coordinates from each numbat run
cnv <- lapply(samples, function(sample.name){
  message(sample.name)
  cnv <- readRDS(str_glue("{data.dir}/{sample.name}/numbat/numbat.rds"))$joint_post %>%
    dplyr::select(CHROM, cnv_state, seg_start, seg_end, seg, LLR) %>%
    distinct() %>%
    mutate(CHROM = str_glue('chr{CHROM}'),
           sample = sample.name)
  return(cnv)
  }) %>%
  do.call(rbind,.)
  #Works now

# Plot --------------------------------------------------------------------

# We can create a plot of bed-like data with the karyoploteR R package's
# 'kpPlotRegions' function. See this link for details:
# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotRegions/PlotRegions.html

# Remove CNVs with a log-likelihood ratio less than 500
#cnv <- filter(cnv, LLR > 500) leaving this out for now.

# Change plotting parameters as needed
pp <- getDefaultPlotParams(plot.type = 3)
pp$leftmargin <- 0.15

# Define some genes of interest and get their genomic coordinates using biomaRt
gene.symbols <- c("CCND1", "FGF4", "CDKN2A", "FGF19", "CDKN2B", "PIK3CA", "MYC")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

#plot
pdf(str_glue("{fig_dir}/allpatients_karyoploteR.pdf"), width=11.6, height=8.2)
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
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)



# Close the plot
dev.off()

print("chromosome 1")


pdf(str_glue("{fig_dir}/allpatients_karyoploteR_Chr1.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr1'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr2")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr2.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr2'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr3")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr3.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr3'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr4")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr4.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr4'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr5")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr5.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr5'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr6")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr6.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr6'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr7")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr7.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr7'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr8")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr8.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr8'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr9")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr9.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr9'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr10")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr10.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr10'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr11")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr11.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr11'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr12")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr12.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr12'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr13")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr13.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr13'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr14")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr14.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr14'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr15")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr15.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr15'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr16")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr16.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr16'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr17")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr17.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr17'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr18")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr18.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr18'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr19")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr19.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr19'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr20")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr20.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr20'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr21")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr21.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr21'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest


# Close the plot
dev.off()

print("chr22")
pdf(str_glue("{fig_dir}/allpatients_karyoploteR_chr22.pdf"), width=11.6, height=8.2)
# Start with an empty plot
kp <- plotKaryotype(genome = "hg38",
                    chromosomes = str_glue('chr22'),
                    plot.type = 3,
                    plot.params = pp)

# Plot amplifications
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "amp")),
              data.panel = 1,
              col = "#cb181d",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)
dev.off()

#Plot each karyoplot with the common genes affected in expanding clones
#Load the expanding clones list:
c<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/common_expanding_genes.csv")
#Create a vector
c<-c$c
#Put into gene.symbols
gene.symbols<-c
#Start the karyoploteR workflow
#Start with HN200519A
sample.name<-"HN200519A"
sample.name<-"HN120819A"
sample.name<-"HN021219A"
sample.name<-"HN230620A"
cnv<-read_csv(file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/{sample.name}_CNV_summary.csv"))

#Define gene symbols
# Define some genes of interest and get their genomic coordinates using biomaRt
gene.symbols <- c
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=108)
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"
#Create the plot
pdf(str_glue("{fig_dir}/{sample.name}_karyoploteR_allchrom_commongenes.pdf"), width=11.6, height=8.2)
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
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="amplification", data.panel = 1)

# Plot deletions
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "del")),
              data.panel = 2,
              col = "#084594",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels="deletion", data.panel = 2, r1 = 0.5)

# Plot LOH
kpPlotRegions(kp,
              data = makeGRangesFromDataFrame(filter(cnv, cnv_state == "loh")),
              data.panel = 2,
              col = "#006d2c",
              r0 = 0.5, r1 = 0.75)
kpAddLabels(kp, labels="LOH", data.panel = 2, r1 = 1)

# Add labels for our genes of interest
kpPlotMarkers(kp,
              data=genes,
              labels=genes$hgnc_symbol,
              line.color = "#555555",
              marker.parts = c(0.95,0.025,0.025),
              r1=0.7)

# Close the plot
dev.off()
#it seems that chromosome 8 is the most likely to be affected.
