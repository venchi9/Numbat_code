#Taking the DEG lists between expanding and non expanding clones and seeing if there are upregulated pathways.
#The data is here:/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists

library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ReactomePA)
library(org.Hs.eg.db)
library(glue)
library(clusterProfiler)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"

#Load in data
sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"
sample_name<-"combined"

x<-read.csv(file = str_glue("{data_dir}/{sample_name}_DEG_list.csv"), row.names =1)

#Insert rownames (SYMBOLS) into a column in the data.frame
x$SYMBOL<-row.names(x)
#Change the gene names to entrez ID
an<-bitr(rownames(x), fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db")
#Match
y<-merge(x, an, by = "SYMBOL")

#I think I will need to filter by p-value
z<-subset(y, y$p_val_adj <= 0.05)

#Select only the entrezid and the log2FC
df<-select(z, c("ENTREZID", "avg_log2FC"))

#Get rid of colnames
colnames(df)<-NULL

#Get rid of rows with NA values
df<-na.omit(df)

#Save
write.table(df, file = str_glue("{data_dir}/{sample_name}_DEG_list_filtered.txt"))

#Worked so far

#Start reactome workflow
geneList = df[,2]
names(geneList) = as.character(df[,1])
geneList = sort(geneList, decreasing = TRUE)
library(ReactomePA)
de <- names(geneList)[abs(geneList) > 0.1]
head(de)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

#savetable
write.csv(x, file = str_glue("{data_dir}/{sample_name}_reactome_enriched_pathways.csv"))

#Now do plots
pdf(str_glue("{fig_dir}/{sample_name}_barplot.pdf"), width=8.2, height=11.6)
barplot(x, showCategory=8)
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_dotplot.pdf"), width=8.2, height=11.6)
dotplot(x, showCategory=15)
dev.off()

pdf(str_glue("{fig_dir}/{sample_name}_cnetplot.pdf"), width=8.2, height=11.6)
cnetplot(x, categorySize="pvalue", foldChange=geneList)
dev.off()

#Worked

#Try merging altogether



