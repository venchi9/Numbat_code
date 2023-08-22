#I have all the correlation matrices now.  I am on a mission to stitch them together.

suppressPackageStartupMessages({
library(SPOTlight)
library(tidyverse)
library(glue)
library(gtools)
library(ggplot2)
})

#Data is here:
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Visium/Seurat/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Visium/Seurat/figures"


print("load in samples")
pt1<-"HN200519A"
pt2<-"HN120819A"
pt3<-"HN021219A"
pt4<-"HN230620A"

print("start with tcells")
group_name <- "tcell"
site<-"primary"

print("Load in matrices")
pt1.cor<-read.csv(file = str_glue("{data_dir}/{pt1}_{group_name}_{site}_cor_matrix.csv"))
pt2.cor<-read.csv(file = str_glue("{data_dir}/{pt2}_{group_name}_{site}_cor_matrix.csv"))
pt3.cor<-read.csv(file = str_glue("{data_dir}/{pt3}_{group_name}_{site}_cor_matrix.csv"))
pt4.cor<-read.csv(file = str_glue("{data_dir}/{pt4}_{group_name}_{site}_cor_matrix.csv"))

print("add meta-data column for each patient")
pt1.cor$patientID<-pt1
pt2.cor$patientID<-pt2
pt3.cor$patientID<-pt3
pt4.cor$patientID<-pt4

print('bind together')
final<-smartbind(pt2.cor, pt3.cor, pt4.cor)

print("extract rows for Expanding and Non.Expanding")
ex<-"Expanding"
#Now try and subet the matrix
ex.table<-final[final$X %in% ex, ]

print("save out")
write.csv(ex.table, file = str_glue("{data_dir}/{site}_{group_name}_{ex}_combined_cor_matrix.csv"))

nonex<-"Non.Expanding"
nonex.table<-final[final$X %in% nonex, ]
write.csv(nonex.table, file = str_glue("{data_dir}/{site}_{group_name}_{nonex}_combined_cor_matrix.csv"))

print("start with bcells")
group_name <- "bcell"
site<-"primary"

print("Load in matrices")
pt1.cor<-read.csv(file = str_glue("{data_dir}/{pt1}_{group_name}_{site}_cor_matrix.csv"))
pt2.cor<-read.csv(file = str_glue("{data_dir}/{pt2}_{group_name}_{site}_cor_matrix.csv"))
pt3.cor<-read.csv(file = str_glue("{data_dir}/{pt3}_{group_name}_{site}_cor_matrix.csv"))
pt4.cor<-read.csv(file = str_glue("{data_dir}/{pt4}_{group_name}_{site}_cor_matrix.csv"))

print("add meta-data column for each patient")
pt1.cor$patientID<-pt1
pt2.cor$patientID<-pt2
pt3.cor$patientID<-pt3
pt4.cor$patientID<-pt4

print('bind together')
final<-smartbind( pt2.cor, pt3.cor, pt4.cor)

print("extract rows for Expanding and Non.Expanding")
ex<-"Expanding"
#Now try and subet the matrix
ex.table<-final[final$X %in% ex, ]
ex.table$clones<-"Expanding"

print("save out")
write.csv(ex.table, file = str_glue("{data_dir}/{site}_{group_name}_{ex}_combined_cor_matrix.csv"))

#Non-expanding
nonex<-"Non.Expanding"
nonex.table<-final[final$X %in% nonex, ]
write.csv(nonex.table, file = str_glue("{data_dir}/{site}_{group_name}_{nonex}_combined_cor_matrix.csv"))

#nonex.table$clones<-"NonExpanding"
#nonex.means<-colMeans(nonex.table[,!colnames(nonex.table) %in% c("X","patientID")], na.rm = T)
#nonex.means<-data.frame(nonex.means)
#nonex.means<-data.frame(t(nonex.means))
#print("Add in meta-data column")
#nonex.means$clones<-"Non-Expanding"

#print("bind expanding and non-expanding means together")
#means<-smartbind(ex.means, nonex.means)

#print("get rid of the Expanding and Non.Expanding meta-data columns")
#means$Expanding<-NULL 
#means$Non.Expanding<-NULL 
3print("might need to transpose for the graphing")
#tmeans<-data.frame(t(means))
#colnames(tmeans)<-c("Expanding", "Non-Expanding")
#remove<-"clones"
#tmeans<-tmeans[!(row.names(tmeans) %in% remove), ]
#tmeans$Celltype<-rownames(tmeans)

#This is easier if we use  a pivot table:
#dfm <- pivot_longer(tmeans, -Celltype, names_to="variable", values_to="value")

#Try to plot
#pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot.pdf"), width=11.6, height=8.2)
#p<-ggplot(dfm, aes(x = Celltype, y = value)) + geom_bar(aes(fill = variable), stat = "identity", position = "dodge") + labs(x = "Celltype", y = "Correlation Interaction") + theme_bw() + theme(axis.text.x = element_text(angle = 45)) + stat_compare_means(aes(group = Celltype), label = "p.signif")
  
#p
#dev.off()
#works.
#I need to work out how to do this with the replicates.
#combined<-smartbind(ex.table, nonex.table)
#combined$patientID<-NULL
#combined$Non.Expanding<-NULL
#combined$Expanding<-NULL
#combined$X<-NULL

#Create pivot table
#                        pivot_longer(!c(clones), names_to = "Celltype", values_to = "value")

###pdf(str_glue("{fig_dir}/{site}_{group_name}_means_bar_plot_pvalues.pdf"), width=11.6, height=8.2)
#p<-ggbarplot(test, x = "Celltype", y = "Correlation", add = "mean_se",
#          color = "clones", palette = "jco", fill = "clones",
#          position = position_dodge(0.8))+
#  stat_compare_means(aes(group = clones), label = "p.signif") + theme(axis.text.x = element_text(angle = 45))
#p
#dev.off()


#Save the combined table
write.csv(combined, str_glue("{data_dir}/{site}_{group_name}_combined_values.csv"))


dfm <- pivot_longer(test, -Celltypes, names_to="variable", values_to="value")

#Save the means table
write.csv(means, file = str_glue("{data_dir}/{site}_{group_name}_combined_means.csv"))



write.csv(nonex.table, file = str_glue("{data_dir}/{site}_{group_name}_{nonex}_combined_cor_matrix.csv"))


print("start with mono")
group_name <- "mono"
site<-"primary"

print("Load in matrices")
pt1.cor<-read.csv(file = str_glue("{data_dir}/{pt1}_{group_name}_{site}_cor_matrix.csv"))
pt2.cor<-read.csv(file = str_glue("{data_dir}/{pt2}_{group_name}_{site}_cor_matrix.csv"))
pt3.cor<-read.csv(file = str_glue("{data_dir}/{pt3}_{group_name}_{site}_cor_matrix.csv"))
pt4.cor<-read.csv(file = str_glue("{data_dir}/{pt4}_{group_name}_{site}_cor_matrix.csv"))

print("add meta-data column for each patient")
pt1.cor$patientID<-pt1
pt2.cor$patientID<-pt2
pt3.cor$patientID<-pt3
pt4.cor$patientID<-pt4

print('bind together')
final<-smartbind( pt2.cor, pt3.cor, pt4.cor)

print("extract rows for Expanding and Non.Expanding")
ex<-"Expanding"
#Now try and subet the matrix
ex.table<-final[final$X %in% ex, ]

print("save out")
write.csv(ex.table, file = str_glue("{data_dir}/{site}_{group_name}_{ex}_combined_cor_matrix.csv"))

nonex<-"Non.Expanding"
nonex.table<-final[final$X %in% nonex, ]
write.csv(nonex.table, file = str_glue("{data_dir}/{site}_{group_name}_{nonex}_combined_cor_matrix.csv"))

print("NK")
group_name <- "NK"
site<-"primary"

print("Load in matrices")
pt1.cor<-read.csv(file = str_glue("{data_dir}/{pt1}_{group_name}_{site}_cor_matrix.csv"))
pt2.cor<-read.csv(file = str_glue("{data_dir}/{pt2}_{group_name}_{site}_cor_matrix.csv"))
pt3.cor<-read.csv(file = str_glue("{data_dir}/{pt3}_{group_name}_{site}_cor_matrix.csv"))
pt4.cor<-read.csv(file = str_glue("{data_dir}/{pt4}_{group_name}_{site}_cor_matrix.csv"))

print("add meta-data column for each patient")
pt1.cor$patientID<-pt1
pt2.cor$patientID<-pt2
pt3.cor$patientID<-pt3
pt4.cor$patientID<-pt4

print('bind together')
final<-smartbind( pt2.cor, pt3.cor, pt4.cor)

print("extract rows for Expanding and Non.Expanding")
ex<-"Expanding"
#Now try and subet the matrix
ex.table<-final[final$X %in% ex, ]

print("save out")
write.csv(ex.table, file = str_glue("{data_dir}/{site}_{group_name}_{ex}_combined_cor_matrix.csv"))

nonex<-"Non.Expanding"
nonex.table<-final[final$X %in% nonex, ]
write.csv(nonex.table, file = str_glue("{data_dir}/{site}_{group_name}_{nonex}_combined_cor_matrix.csv"))

print("other")
group_name <- "other"
site<-"primary"

print("Load in matrices")
pt1.cor<-read.csv(file = str_glue("{data_dir}/{pt1}_{group_name}_{site}_cor_matrix.csv"))
pt2.cor<-read.csv(file = str_glue("{data_dir}/{pt2}_{group_name}_{site}_cor_matrix.csv"))
pt3.cor<-read.csv(file = str_glue("{data_dir}/{pt3}_{group_name}_{site}_cor_matrix.csv"))
pt4.cor<-read.csv(file = str_glue("{data_dir}/{pt4}_{group_name}_{site}_cor_matrix.csv"))

print("add meta-data column for each patient")
pt1.cor$patientID<-pt1
pt2.cor$patientID<-pt2
pt3.cor$patientID<-pt3
pt4.cor$patientID<-pt4

print('bind together')
final<-smartbind( pt2.cor, pt3.cor, pt4.cor)

print("extract rows for Expanding and Non.Expanding")
ex<-"Expanding"
#Now try and subet the matrix
ex.table<-final[final$X %in% ex, ]

print("save out")
write.csv(ex.table, file = str_glue("{data_dir}/{site}_{group_name}_{ex}_combined_cor_matrix.csv"))

nonex<-"Non.Expanding"
nonex.table<-final[final$X %in% nonex, ]
write.csv(nonex.table, file = str_glue("{data_dir}/{site}_{group_name}_{nonex}_combined_cor_matrix.csv"))

