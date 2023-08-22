#Find an intersected gene list from each enriched pathway for the DEG analysis with clones


fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(cowplot)
library(tidyverse)
library(glue)

#Data is here:
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
list_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/reactomeGeneLists"

#The commonly enriched pathways are:

#Viral mRNA Translation
#Signaling by ROBO receptors
#Regulation of expression of SLITs and ROBOs
#Formation of a pool of free 40S subunits
#Peptide chain elongation

#Samples
sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"

#Comparisons
comparison_name<-"LN_v_primary4"
comparison_name<-"LN_v_primary6" - no viral mRNA translation
comparison_name<-"LN2_vP2"
comparison_name<-"LN3_v_P3"
comparison_name<-"LN23_v_P23"
comparison_name<-"LN3_v_primary"
comparison_name<-"LN4_v_primary"
comparison_name<-"LN_v_primary"

#Lists are here:
#file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#a = HN200519A LNv primary4
#b = HN200519A LN vprimary 6
#c = HN120819A LN2 vP2
#d = HN120819A LN3 v P3
#e = HN120819A LN23 vP23
#f = HN021219A LN3 v primary
#g = HN021219A LN4 v primary
#h = HN230620A LN v primary

#Start with a
sample_name<-"HN200519A"
comparison_name<-"LN_v_primary4"
a<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#b
sample_name<-"HN200519A"
comparison_name<-"LN_v_primary6"
b<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#c
sample_name<-"HN120819A"
comparison_name<-"LN2_vP2"
c<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#d
sample_name<-"HN120819A"
comparison_name<-"LN3_v_P3"
d<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))
#e
sample_name<-"HN120819A"
comparison_name<-"LN23_v_P23"
e<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#f
sample_name<-"HN021219A"
comparison_name<-"LN3_v_primary"
f<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#g
sample_name<-"HN021219A"
comparison_name<-"LN4_v_primary"
g<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#h
sample_name<-"HN230620A"
comparison_name<-"LN_v_primary"
h<-read.csv(file = str_glue("{list_dir}/{sample_name}_{comparison_name}_reactome_pathway_list.csv"))

#Look for intesecting pathways
#Viral mRNA Translation
#Signaling by ROBO receptors
#Regulation of expression of SLITs and ROBOs
#Formation of a pool of free 40S subunits
#Peptide chain elongation

#Start with Viral.mRNA.Translation
#Create a vector for each comparison
#a
aviral<-a$Viral.mRNA.Translation
#b
nul
#c
cviral<-c$Viral.mRNA.Translation
#d
dviral<-d$Viral.mRNA.Translation
#e
eviral<-e$Viral.mRNA.Translation
#f
fviral<-f$Viral.mRNA.Translation
#g
Nil
#h
hviral<-h$Viral.mRNA.Translation

#So have lists in a, c, d, e, f, h
virallist<-Reduce(intersect, list(aviral, cviral, dviral, eviral, fviral, hviral))
write.csv(virallist, file = str_glue("{list_dir}/viralmRNA_intersection_list.csv"))

#Signaling.by.ROBO.receptors
#a
arobo<-a$Signaling.by.ROBO.receptors
#b
brobo<-b$Signaling.by.ROBO.receptors
#c
crobo<-c$Signaling.by.ROBO.receptors
#d
drobo<-d$Signaling.by.ROBO.receptors
#e
erobo<-e$Signaling.by.ROBO.receptors
#f
nul
#g
grobo<-g$Signaling.by.ROBO.receptors
#h
Nul

#So have lists in a,b,c,d,e,g
robolist<-Reduce(intersect, list(arobo, brobo, crobo, drobo, erobo, grobo))
write.csv(virallist, file = str_glue("{list_dir}/ROBO_intersection_list.csv"))

#Regulation.of.expression.of.SLITs.and.ROBOs
#a
aslit<-a$Regulation.of.expression.of.SLITs.and.ROBOs
#b
bslit<-b$Regulation.of.expression.of.SLITs.and.ROBOs
#c
cslit<-c$Regulation.of.expression.of.SLITs.and.ROBOs
#d
dslit<-d$Regulation.of.expression.of.SLITs.and.ROBOs
#e
eslit<-e$Regulation.of.expression.of.SLITs.and.ROBOs
#f
nul
#g
gslit<-g$Regulation.of.expression.of.SLITs.and.ROBOs
#h
nul

slitlist<-Reduce(intersect, list(aslit, bslit, cslit,dslit, eslit, gslit))
write.csv(slitlist, file = str_glue("{list_dir}/SLIT_intersection_list.csv"))


#Formation.of.a.pool.of.free.40S.subunits

#a
asub<-a$Formation.of.a.pool.of.free.40S.subunits
#b
bsub<-b$Formation.of.a.pool.of.free.40S.subunits
#c
csub<-c$Formation.of.a.pool.of.free.40S.subunits
#d
dsub<-d$Formation.of.a.pool.of.free.40S.subunits
#e
esub<-e$Formation.of.a.pool.of.free.40S.subunits
#f
fsub<-f$Formation.of.a.pool.of.free.40S.subunits
#g
nul
#h
hsub<-h$Formation.of.a.pool.of.free.40S.subunits


sublist<-Reduce(intersect, list(asub, bsub, csub, dsub, esub, fsub, hsub))
write.csv(sublist, file = str_glue("{list_dir}/40S_intersection_list.csv"))

#Peptide.chain.elongation
apep<-a$Peptide.chain.elongation
#b
bpep<-b$Peptide.chain.elongation
#c
cpep<-c$Peptide.chain.elongation
#d
dpep<-d$Peptide.chain.elongation
#e
epep<-e$Peptide.chain.elongation
#f
fpep<-f$Peptide.chain.elongation
#g
nul
#h
hpep<-h$Peptide.chain.elongation
peplist<-Reduce(intersect, list(apep, bpep, cpep, dpep, epep, fpep, hpep))
write.csv(peplist, file = str_glue("{list_dir}/Peptide_intersection_list.csv"))


