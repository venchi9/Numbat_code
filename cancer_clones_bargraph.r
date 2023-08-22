#Take the cancer clones and make bargraphs for figure 2 of the paper
fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(glue)
library(gtools)
library(ggplot2)
library(ggpubr)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/figures"

sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"


print("load in the data")

clones<-read.csv(file = str_glue("{data_dir}/{sample_name}_clone_percentages.csv"), header=TRUE)

print("create a pivot table")

require(reshape2)
melted<-melt(data = clones, id.vars = 1)
names(melted)<-c("Clones", "Site", "Percentage")
melted$Clones<-as.factor(melted$Clones)
print(melted)

print("try plotting now")
pdf(str_glue("{fig_dir}/{sample_name}_clone_percentages_by_clone.pdf"), width=8.2, height=11.6)
p<-ggplot(data = melted, aes(x = Clones, y = Percentage, fill = Site)) + geom_bar(stat = "identity", position = position_dodge())
p + theme_bw()+ ylim(ylim = c(0, 100))
dev.off()



print("plot bars by site")
pdf(str_glue("{fig_dir}/{sample_name}_clone_percentages_by_site.pdf"), width=8.2, height=11.6)
p<-ggplot(data = melted, aes(x = Site, y = Percentage, fill = Clones)) + geom_bar(stat = "identity", position = position_dodge())
p + scale_fill_manual(values = c("#BEBEBE", "#E41A1C", "#3A86A5", "#658E67", "#CB6651", "#FFD422", "#B47229")) + theme_bw() + ylim(ylim = c(0, 100))
dev.off()

print("plot stacked bar graph by site")
pdf(str_glue("{fig_dir}/{sample_name}_clone_stackedbar_bysite.pdf"), width=8.2, height=11.6)
p<-ggplot(data = melted, aes(x = Site, y = Percentage, fill = Clones)) + geom_bar(stat="identity")
p + scale_fill_manual(values = c("#BEBEBE", "#E41A1C", "#3A86A5", "#658E67", "#CB6651", "#FFD422", "#B47229"))+ theme_bw() + ylim(ylim = c(0, 100))
dev.off()


