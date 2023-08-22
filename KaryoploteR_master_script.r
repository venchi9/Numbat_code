#Do KaryoploteR on all the samples for the paper.

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

pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_allchromosomes.pdf"), width=11.6, height=8.2)
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


# Close the plot
dev.off()

print("chromosome 1")


pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_Chr1.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr2.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr3.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr4.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr5.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr6.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr7.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr8.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr9.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr10.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr11.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr12.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr13.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr14.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr15.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr16.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr17.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr18.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr19.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr20.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr21.pdf"), width=11.6, height=8.2)
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
pdf(str_glue("{fig_dir}/{sample_name}_karyoploteR_chr22.pdf"), width=11.6, height=8.2)
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

# Add labels for our genes of interest


# Close the plot
