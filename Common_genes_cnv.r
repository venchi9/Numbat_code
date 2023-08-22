#Look for common genes affected by CNV 
#Data is here:
/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate r_4

library(org.Hs.eg.db)
library(GenomicFeatures)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(Homo.sapiens)

data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures"
code_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/code"

#Make a geneRanges function
geneRanges <- 
    function(db, column="ENTREZID")
{
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
}

splitColumnByOverlap <-
    function(query, subject, column="ENTREZID", ...)
{
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}

#Load in data:
cnv<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/combined_cnv.csv")
#Tidy up a bit
cnv$...1<-NULL
cnv$...2<-NULL
colnames(cnv)[3]<-"start"
colnames(cnv)[4]<-"end"

#Make into GRanges output
cnv = makeGRangesFromDataFrame(cnv)

library(Homo.sapiens)

gns = geneRanges(Homo.sapiens, column="SYMBOL")
#Find which genes overlap in which regions
symInCnv = splitColumnByOverlap(gns, cnv, "SYMBOL")
symInCnv
list<-as.data.frame(symInCnv)
#Clean up
list$group<-NULL
list$group_name<-NULL
colnames(list)[1]<-"gene"

#save this list

write.csv(list, file = str_glue("{data_dir}/genelist_fromCNV.csv"))
#Create vector with genes listed

******************************
#Himanshi did this for me and the data is here:
#/directflow/SCCGGroupShare/projects/himaro/projects/head_neck_venessa

#Read the genes which had an ensemble ID and were seen in the orginal sc matrix
x<-read_tsv(file = "/directflow/SCCGGroupShare/projects/himaro/projects/head_neck_venessa/new_dt_matched_genes.tsv")
#dim = 32114
#Read in genes which did not match the sc matrix
y<-read_tsv(file = "/directflow/SCCGGroupShare/projects/himaro/projects/head_neck_venessa/mismatch_genes.tsv")
#dim = 5360 - looks like a lot of pseudogenes

#Create a list of the object x of the ensembl gene ids
ensembl_gene_ids<-x$ensebl_id

#Load up biomaRt
library(biomaRt)

listEnsemblArchives()

ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl",
                      version = 108)

#Create a list of the ensembl gene ids and their classification
list <- getBM(attributes = c("gene_biotype", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = ensembl_gene_ids,
      mart = ensembl)

# head(x)
                            lncRNA                     protein_coding 
                               987                              14897 
  transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
                                 1                                  5 
transcribed_unprocessed_pseudogene             unprocessed_pseudogene 
                                11                                  1 

#Subset the list to those which are protein_coding

x<-subset(x, subset = x$gene == "protein_coding")
list<-subset(list, subset = list$gene_biotype == "protein_coding")

#Okay, I think I'll do this by patient and can then do a VENN diagram to see the overlap.
#HN200519A

#Load in data:
cnv<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/HN200519A_CNV_summary.csv")
cnv<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/HN120819A_CNV_summary.csv")
cnv<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/HN021219A_CNV_summary.csv")
cnv<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/HN230620A_CNV_summary.csv")
#Tidy up a bit
cnv$...1<-NULL
cnv$...2<-NULL
colnames(cnv)[3]<-"start"
colnames(cnv)[4]<-"end"

#Make into GRanges output
cnv = makeGRangesFromDataFrame(cnv)

library(Homo.sapiens)

gns = geneRanges(Homo.sapiens, column="SYMBOL")
#Find which genes overlap in which regions
symInCnv = splitColumnByOverlap(gns, cnv, "SYMBOL")
symInCnv
list<-as.data.frame(symInCnv)
#Clean up
list$group<-NULL
list$group_name<-NULL
colnames(list)[1]<-"gene"

#Create a vector list of the genes in HN200519A
HN200519A<-list
HN200519A<-HN200519A$gene
#Create vector list of genes in HN120819A
HN120819A<-list
HN120819A<-HN120819A$gene

#Create vector list of genes in HN021219A
HN021219A<-list
HN021219A<-HN021219A$gene

#Create vector list of genes in HN230620A
HN230620A<-list
HN230620A<-HN230620A$gene

#Now do the intersection
a<-intersect(HN200519A, HN120819A)
b<-intersect(HN021219A, HN230620A)
c<-intersect(a, b)
#length is 378

#Load in the data with the ensembl ids
x<-read_tsv(file = "/directflow/SCCGGroupShare/projects/himaro/projects/head_neck_venessa/new_dt_matched_genes.tsv")
#Subset this to the genes in c

final<-subset(x, subset = x$gene %in% c)
#Return only unique values
test<-distinct(final, gene, .keep_all = TRUE)
test<-na.omit(na)
#Final number is 317

#Create a list of the object x of the ensembl gene ids
ensembl_gene_ids<-test$ensebl_id

#Load up biomaRt
library(biomaRt)

listEnsemblArchives()

ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl",
                      version = 108)

#Create a list of the ensembl gene ids and their classification
list <- getBM(attributes = c("gene_biotype", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = ensembl_gene_ids,
      mart = ensembl)

# head(x)
                            lncRNA                     protein_coding 
                               987                              14897 
  transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
                                 1                                  5 
transcribed_unprocessed_pseudogene             unprocessed_pseudogene 
                                11                                  1 

#Subset the list to those which are protein_coding

finallist<-subset(list, subset = list$gene_biotype == "protein_coding")
#Get list of genes
genes<-finallist$ensembl_gene_id

#Subset original dataframe for the protein coding genes
v<-subset(test, subset = test$ensebl_id %in% genes)
#dim = 281

#save out
write.csv(v, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/commonly_affected_coding_genes.csv")


#Next create master list of all the protein coding genes from combined list and then subset the genes from each patient --> make a VENN diagram
#Also see if these genes have any commonalities, use stringDB
#First load up the master list
x<-read_tsv(file = "/directflow/SCCGGroupShare/projects/himaro/projects/head_neck_venessa/new_dt_matched_genes.tsv")
#Then create the list of ensembl ids which are protein coding genes
ensembl_list<-list$ensembl_gene_id
#Subset the master list to protein coding only
x<-subset(x, subset = x$ensebl_id %in% ensembl_list)
#Keep only distinct values
x<-distinct(x, ensebl_id, .keep_all = TRUE)
#1496
#Save out
write.csv(x, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/karyoploter/data/protein_coding_genes_combined.csv")
#Make a vector
master<-x$gene
#Load in individual lists
#Now do the intersections with the master list which will just give us the coding genes

a<-intersect(HN200519A, master)
#5956
b<-intersect(HN120819A, master)
#4449
c<-intersect(HN021219A, master)
#8691
d<-intersect(HN230620A, master)
#5715

#make all vectors the same length
#Longest is 8691

n<-8691
length(HN200519A)<-n
length(HN230620A)<-n
length(HN120819A)<-n
length(HN021219A)<-n

#Try VENN diagram
library("ggVennDiagram")
listInput<-list(patient_1 =  HN200519A, patient_2 = HN120819A, patient_3 = HN021219A,   patient_4 =  HN230620A)
pdf(str_glue("{fig_dir}/common_cnv_genes_VENN_diagram.pdf"), width=11.6, height=8.2)
ggVennDiagram(listInput)
dev.off()
#Not working
library("ggvenn")
  pdf(str_glue("{fig_dir}/common_cnv_genes_VENN_diagram.pdf"), width=11.6, height=8.2)
  listInput<-list(patient_1 =  HN200519A, patient_2 = HN120819A, patient_3 = HN021219A,   patient_4 =  HN230620A)
  ggvenn(listInput)
  dev.off()








