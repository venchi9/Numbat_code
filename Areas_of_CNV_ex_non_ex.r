#Try and get genes from areas of CNV in expanding and non-expanding clones
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

#Samples
  sample.name<-"HN021219A"
  sample.name<-"HN200519A"
  sample.name<-"HN120819A"
  sample.name<-"HN230620A"

#nb object generation:
nb = Numbat$new(str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}/numbat"))

#Of course, HN200519A wont work properly despite my best efforts.  We are using numbat_400 for the clones
x1<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN200519A/numbat_400/joint_post_20.tsv")
x2<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN200519A/numbat_400/clone_post_20.tsv")


#Of course, Walter cracked this.

x1 <- dplyr::select(nb$joint_post, cell, CHROM, seg, cnv_state,seg_start, seg_end)
x2 <- dplyr::select(nb$clone_post, cell, clone_opt)
x3 <- merge(x1, x2, by = "cell")
x3 <- dplyr::select(x3, -cell)
x3 <- distinct(x3)


#In HN021219A, clones 3 and 4 are expanding
x3$exclones<-ifelse((x3$clone_opt == "3") | (x3$clone_opt == "4"), "expanding", "other")

#IN HN200519A, clone 3 and 5 are expanding
x3$exclones<-ifelse(((x3$clone_opt == "3") | (x3$clone_opt == "5")), "expanding", "other")
#IN HN120819A, clone 3 is expanding
x3$exclones<-ifelse((x3$clone_opt == "3"), "expanding", "other")
#IN HN230620A, clones 4 and 5 are expanding
x3$exclones<-ifelse((x3$clone_opt == "4") | (x3$clone_opt == "5"), "expanding", "other")

#Separate out the ex and nonex tables for saving.
expanding<-subset(x3, subset = x3$exclones =="expanding")
#Because there are lots of shared CNVs with the clones, we need to find the distinct ones.
#Remove clone informatio

expanding<-dplyr::select(expanding, -clone_opt)

#Get unique values
expanding<-distinct(expanding)

#save out
write.csv(expanding, file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}_expanding_list.csv"))
write.csv(nonex, file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}_nonexpanding_list.csv"))

#Now need to work on this data more.  Firstly, create a list of genes affected by CNV for each patient, ex and non-ex.
library(org.Hs.eg.db)
library(GenomicFeatures)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(Homo.sapiens)

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

#Read in masterlist of coding genes
#Read the genes which had an ensemble ID and were seen in the orginal sc matrix
x<-read_tsv(file = "/directflow/SCCGGroupShare/projects/himaro/projects/head_neck_venessa/new_dt_matched_genes.tsv")
#Find from ids which are coding
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

#Change table 'x' colnames
colnames(x)<-c("gene", "ensemble_gene_ids")
#Merge tables together
list<-subset(list, subset = list$gene_biotype == "protein_coding")
coding<-list$ensembl_gene_id
#Subset x by vector coding
master<-subset(x, subset = x$ensemble_gene_ids %in% coding)
#remove NA values
master<-na.omit(master)
#Distinct values only
master<-distinct(master, ensemble_gene_ids, .keep_all = TRUE)
#Save the master list of coding genes
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones"
write.csv(master, file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/master_coding_gene_list.csv"))

#Start with patient HN200519A
#Load in data.

  sample.name<-"HN021219A"
  sample.name<-"HN200519A"
  sample.name<-"HN120819A"
  sample.name<-"HN230620A"

cnv<-read.csv(file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}_expanding_list.csv"))
cnv$X<-NULL
cnv$CHROM<-paste0("chr", cnv$CHROM)

#Change column names
colnames(cnv)[4]<-"start"
colnames(cnv)[5]<-"end"


#Keep going for now and see what happens

#Make into GRanges output
cnv = makeGRangesFromDataFrame(cnv)

library(Homo.sapiens)

gns = geneRanges(Homo.sapiens, column="SYMBOL")
#Find which genes overlap in which regions
symInCnv = splitColumnByOverlap(gns, cnv, "SYMBOL")
symInCnv
cnvlist<-as.data.frame(symInCnv)
#Clean up
cnvlist$group<-NULL
cnvlist$group_name<-NULL
colnames(cnvlist)[1]<-"gene"

#Create a vector of genes
cnvvec<-cnvlist$gene

#

#Filter for coding genes
mastervec<-master$gene

#Filter expanding for coding genes
cnvcoding<-intersect(mastervec, cnvvec)

cnvcoding<-as.data.frame(cnvcoding)

#Saveout
write.csv(cnvcoding, file = str_glue("/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/{sample.name}_expanding_coding_genes.csv"))


#Now load up all the different patients
HN2005<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN200519A_expanding_coding_genes.csv")
#5956
HN0212<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN021219A_expanding_coding_genes.csv")
#8691
HN1208<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN120819A_expanding_coding_genes.csv")
#4449
HN2306<-read_csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/HN230620A_expanding_coding_genes.csv")
#5715

#Now find the common genes between them.
#Firstly, create vectors
HN2005<-HN2005$cnvcoding
HN0212<-HN0212$cnvcoding
HN1208<-HN1208$cnvcoding
HN2306<-HN2306$cnvcoding

a<-intersect(HN2005, HN0212)
b<-intersect(HN1208, HN2306)
c<-intersect(a, b)
#281 genes
c<-as.data.frame(c)
#Write this out
write.csv(c, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Numbat/data/common_expanding_genes.csv")