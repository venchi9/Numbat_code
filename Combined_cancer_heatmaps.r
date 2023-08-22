#Heatmaps for the love of god.

#Trying to make sense of the expression data.  I have two concepts here.  The protein translation pathway and the immune pathway.

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(cowplot)
library(tidyverse)
library(glue)
library(Seurat)
library(scCustomize)
library(stats)
library(rstatix)
library(ggpubr)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
list_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/reactomeGeneLists"

allmerged<-readRDS(file = str_glue("{data_dir}/combined_cancer_allmerged_allfeatures.rds"))

#Find All Markers
Idents(allmerged)<-"final"
markers<-FindAllMarkers(allmerged, min.pct = 0.1, logfc.threshold = 0.25)
#Filter by p value
markers1<-filter(markers, p_val_adj <0.05)
#Load in the protein translation list
list<-read.csv(file = str_glue("{data_dir}/combined_protein_trans_list.txt"))
#Create vector
list<-list$x
#Intersect with variable genes
variable<-markers1$gene
#Find intersection
trans_list<-intersect(list, variable)

#Now, limit the markers1 data table by the trans_list
trans.df<-subset(markers1, markers1$gene %in% trans_list)
#Order the list
top20 <- trans.df %>%group_by(cluster) %>%top_n(71, avg_log2FC)
list<-top20$gene
list<-sort(list)
#Do heatmap with clustering
Idents(allmerged)<-"final"
allmerged@active.ident <- factor(allmerged@active.ident, 
                            levels=c("Primary", "LN"))
pdf(str_glue("{fig_dir}/new_combined_cancer_translation_markergenes.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = list)+ NoLegend() + theme(axis.text.y = element_text(size=5.5))
p 
dev.off()

#Not that interesting, but okay I guess

#Now try for the immune genes.
IFN<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/IFN_stimulated_genes.txt", col_names = FALSE)
IFN<-IFN$X1
immune_list<-intersect(IFN, variable)
immune.df<-subset(markers1, markers1$gene %in% immune_list)
top20 <- immune.df %>%group_by(cluster) %>%top_n(9, avg_log2FC)

#Do heatmap with clustering
Idents(allmerged)<-"final"
allmerged@active.ident <- factor(allmerged@active.ident, 
                            levels=c("Primary", "LN"))
pdf(str_glue("{fig_dir}/new_combined_cancer_immune_markergenes.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = top20$gene)+ NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()


#Try an Ag presentation gene list
Ag<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Ag_presentation2.txt", col_names = FALSE)
Ag<-Ag$X1
Ag_list<-intersect(Ag, variable)
ag.df<-subset(markers1, markers1$gene %in% Ag_list)
top20 <- ag.df %>%group_by(cluster) %>%top_n(50, avg_log2FC)
pdf(str_glue("{fig_dir}/new_combined_cancer_Ag_presenation_markergenes.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = top20$gene)+ NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

#Cancer dependent genes
gene_list<-c("ERCC6L", "CKAP2", "CCNA2", "MCM7", "VEGFA", "HIF1A", "SP1", "PCBP1", "PCBP2", "TP53", "BAK1", "FAS", "BAX")
pdf(str_glue("{fig_dir}/new_combined_critical_cancer_genes.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = gene_list)+ NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()
#https://www.cell.com/cell-reports/pdfExtended/S2211-1247(23)00519-3 inducible gene signature

gene_list<-c("OASL", "DDX60", "RSAD2", "IFIT3", "IFIT1", "IFI44", "MX1", "DDX58", "ISG15", "IFIT2", "IFI27", "XAF1", "IFI6", "CMPK2", "MX2", "IFI44L", "OAS3", "HELZ2", "OAS2", "SAMHD1", "HERC6",  "TRIM22", "HERC5", "OAS1", "EPSTI1", "CXCL11", "GBP4", "ZC3HAV1", "DDX60L", "HLA-B", "IRF7", "STAT2", "CXCL10", "DHX58", "SAMD9", "PARP14", "ZNFX1", "AXL")
pdf(str_glue("{fig_dir}/new_combined_IFNresponse_genes.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = gene_list)+ NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()

Idents(allmerged)<-"patients"
Idents(allmerged)<-"patients"
allmerged@active.ident <- factor(allmerged@active.ident, 
                            levels=c("Patient1_Primary",  "Patient2_Primary", "Patient3_Primary",  "Patient4_Primary",  "Patient1_LN", "Patient2_LN","Patient3_LN","Patient4_LN"))

pdf(str_glue("{fig_dir}/new_combined_cancer_translation_markergene_bypatient.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = list)+ NoLegend() + theme(axis.text.y = element_text(size=5.5))
p 
dev.off()

pdf(str_glue("{fig_dir}/new_combined_IFNresponse_genes_bypatient.pdf"), width=11.6, height=8.2)
p<-DoHeatmap(allmerged, features = gene_list)+ NoLegend() + theme(axis.title.y = element_text(size=2))
p 
dev.off()




#Get expression data
gene<-"EIF4G1"

groupA<-subset(allmerged, subset = final =="Primary")
groupB<-subset(allmerged, subset = final =="LN")

expression_A<-data.frame(expression = groupA@assays$SCT@scale.data[gene,])
expression_B<-data.frame(expression = groupB@assays$SCT@scale.data[gene,])

#Combine the data
expression_A$Group<-"Primary"
expression_B$Group<-"LN"
combined<-rbind(expression_A, expression_B)

#Plot

y= combined$expression
x = combined$Group

coefs<-coef(lm(y~x, data = combined))


pdf(str_glue("{fig_dir}/boxplot_{gene}_witht_test.pdf"), width=4, height=4)
level_order<-c('Primary', 'LN')
my_comparisons<-list(c("Primary", "LN"))
p <- ggplot(combined, aes( x = factor(Group, level = level_order), y = expression)) +
  geom_boxplot(aes(fill = Group)) + ylim(-5, 6)
 p + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") + theme_bw()
 dev.off()

 pdf(str_glue("{fig_dir}/boxplot_{gene}_nolegends.pdf"), width=3, height=3)
level_order<-c('Primary', 'LN')
my_comparisons<-list(c("Primary", "LN"))
p <- ggplot(combined, aes( x = factor(Group, level = level_order), y = expression)) +
  geom_boxplot(aes(fill = Group)) 
 p  + theme_bw() + NoLegend() + theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_fill_manual(values = c("red", "darkred"))
 dev.off()

#Conduct a t-test
gene<-"IFIT3"

group1 <- allmerged$final == "Primary"
group2 <- allmerged$final == "LN"

t_test<-t.test(allmerged@assays$SCT@scale.data[gene,group1], allmerged@assays$SCT@scale.data[gene, group2])
t_test