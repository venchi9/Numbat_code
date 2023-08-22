#Take the celltype correlations between the cancer cells and immune cells across the different patients.  Calculate the means and then graph and do p-values for the difference.
#The data is here: /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Visium/Seurat/data

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(SPOTlight)
library(tidyverse)
library(glue)
library(gtools)
library(ggplot2)
library(ggpubr)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Visium/Seurat/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Visium/Seurat/figures"


print("start with tcells")
group_name <- "tcell"
print("start with primary")
site<-"primary"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Get rid of any celltypes with <2 replicates
combined$CD8Naive<-NULL
combined$CD8Proliferating<-NULL
combined$dnT<-NULL 


print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this

pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))

p 
dev.off()

#Works

#Now for the LN
site<-"LN"


print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X<-NULL
combined$X.1<-NULL

#Get rid of celltypes with only 1 replicate
combined$CD8Proliferating<-NULL 

print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#Bcells
site<-"primary"
group_name<-"bcell"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1)) 
p
dev.off()

#LN

site<-"LN"
group_name<-"bcell"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#Monocytes

site<-"primary"
group_name<-"mono"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Remove celltypes with only 1 replicate
combined$pDC<-NULL

print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#LN
site<-"LN"
group_name<-"mono"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Remove celltypes with only 1 replicate
combined$pDC<-NULL

print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#NK
site<-"primary"
group_name<-"NK"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Remove celltypes with only 1 replicate


print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#LN
site<-"LN"
group_name<-"NK"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Remove celltypes with only 1 replicate
combined$NKProliferating<-NULL

print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#Other
site<-"primary"
group_name<-"other"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Remove celltypes with only 1 replicate


print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()

#LN

site<-"LN"
group_name<-"other"

print("load in the data")
ex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Expanding_combined_cor_matrix.csv"))
ex.table$clones<-"Expanding"
nonex.table<-read.csv(file = str_glue("{data_dir}/{site}_{group_name}_Non.Expanding_combined_cor_matrix.csv"))
nonex.table$clones<-"Non-Expanding"
print("combine the tables")
combined<-smartbind(ex.table, nonex.table)
combined$patientID<-NULL
combined$Non.Expanding<-NULL
combined$Expanding<-NULL
combined$X.1<-NULL
combined$X<-NULL

#Remove celltypes with only 1 replicate


print("create pivot table")
#Create pivot table
test <- combined %>% 
                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")


#Try and plot this
pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
p<-ggbarplot(test, x = "Celltype", y = "value", add = "mean_se",
          color = "clones", palette = c("#0073C2","#E46726"), fill = "clones", xlab = FALSE, ylab = FALSE, width = 0.7,
           position = position_dodge(0.8))+
  stat_compare_means(aes(group = clones), label = "p.signif", method = "t.test") + theme(axis.text.x = element_text(angle = 45))+ ylim (ylim = c(-1, 1))
p
dev.off()
