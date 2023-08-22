#Do a summary heatmap for the DEG expression figure 3

fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate r_4
cd /directflow/SCCGGroupShare/projects/venchi/Lung/Never_smokers/eQTL/interactions
R

library(ggplot2)
library(tidyverse)
library(cowplot)
library(glue)
library(ComplexHeatmap)

fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"

#The summary table I made is here:
#/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists/reactome_data/Reactome_summary_sats.csv

test<-read.csv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists/reactome_data/Reactome_summary_sats.csv")
#New data
test<-read.csv(file ="/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/New_transcriptionalgroups_reactome_enriched.csv")

#Change the column names
colnames(test)[2] = "Patient"
colnames(test)[1] = "Pathway"

#Set order of the pathways:
level_order<-c("Collagen Synthesis", "Collagen formation", "ECM degradation", "Collagen degradation", "ECM organisation", "Collagen assembly", "Integrin Interactions","ECM proteoglycans", 
"IL4/IL-13 signaling", "MHC II Ag Presentation", "PD-1 Signaling", "IF-gamma Signaling", "CD3/TCR-zeta phosphorylation", "ZAP-70 to immunological synapse", "IL signaling", 
"Chemokine binding", "IL-10 signaling", "IGF regulation", "Protein phosphorylation", "IGF transport/uptake") 

#New data:
level_order<-c("Translation",
"Peptide chain elongation",
"L13a-mediated translational silencing of Ceruloplasmin expression",
"Formation of a pool of free 40S subunits",
"Eukaryotic Translation Termination",
"Eukaryotic Translation Initiation",
"Eukaryotic Translation Elongation",
"GTP hydrolysis and joining of the 60S ribosomal subunit",
"SRP-dependent cotranslational protein targeting to membrane",
"Cap-dependent Translation Initiation",
"Regulation of expression of SLITs and ROBOs",
"Signaling by ROBO receptors",
"Viral mRNA Translation",
"Response of EIF2AK4 to amino acid deficiency",
"Nonsense Mediated Decay independent of the Exon Junction Complex",
"The citric acid cycle and respiratory electron transport",
"Respiratory electron transport and ATP synthesis",
"Selenocysteine synthesis",
"Interferon Signaling",
"Interferon gamma signaling",
"Diseases of signal transduction by growth factor receptors and second messengers",
"Mitotic Anaphase",
"Mitotic Metaphase and Anaphase",
"Programmed Cell Death",
"Apoptosis",
"Transcriptional Regulation by TP53")

#Open the plot
#Need to pivot the table longer and wider
longdf<-test %>%
    pivot_wider(names_from = Pathway,
                values_from = P.value) %>%
    pivot_longer(
         -Patient,
         names_to = c("Pathway")
       )
pdf(str_glue("{fig_dir}/test_summary_heatmap_newdata.pdf"), width = 8.2, heigh = 11.6)

#Plot

ggplot(longdf, aes(x = Patient, y = factor(Pathway, level = level_order), fill = value)) + geom_tile(color = "grey")+ coord_fixed()+ scale_fill_continuous(na.value = "white") + scale_x_discrete(position = "top")

dev.off()

#It's finally right!





     