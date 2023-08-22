fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(cowplot)
library(tidyverse)
library(glue)

gene_data<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists"
reactome_data<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/genelists/reactome_data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"

sample_name<-"combined"
sample_name<-"HN200519A"
sample_name<-"HN021219A"
sample_name<-"HN120819A"
sample_name<-"HN230620A"

#Step one.  Organise data for gene.set.list:  a named list of the genes in each gene set
#Load up the reactome pathways.
geneset<-read.csv(file = str_glue("{reactome_data}/{sample_name}_reactome_enriched_pathways.csv"))

#Make vectors for the top 10 genelists
x<-geneset[1:10, ]

#Walter's solution
gene.sets <- strsplit(as.character(x$geneID), "/")
names(gene.sets) <- x$Description

#Now get a list of the log2FC for each gene in our data
FC<-read.csv(file = str_glue("{gene_data}/{sample_name}_DEG_list.csv"))
#Filter by significance
z<-subset(FC, FC$p_val_adj <= 0.05)
#Subset out just gene name and log2FC
fold.change<-subset(z, select = c("X", "avg_log2FC"))
#Rename columns
colnames(fold.change)[1] ="gene"
colnames(fold.change)[2] ="log2FC"

#Make the pvals
pvals<-subset(x, select=c("Description", "p.adjust"))
#Change names
colnames(pvals)[1] = "gene_set"
colnames(pvals)[2] = "gsa_pval"

# Define plotting function ------------------------------------------------

# Define a function with the following three arguments as input;
# gene.set.list - a named list of the genes in each gene set
# fold.change.results - a data frame with two columns; `genes` and `log2FC`
# pval.results - a data frame with two columns; `gene_set` and `gsa_pval`
# the function returns a plot with the number of genes up/down regulated for
# each result



gsa_plot <- function(gene.set.list = gene.sets,
                     fold.change.results = fold.change,
                     pval.results = pvals){
    
    # Change inputs into a single data frame for plotting
    plot_data <- map_df(.x = 1:10, .f = function(i) expand.grid(unlist(gene.set.list[[i]]), names(gene.set.list)[i])) %>%
        rename(gene = Var1, gene_set = Var2) %>%
        merge(fold.change.results, by="gene") %>%
        merge(pval.results, by="gene_set")
    
    # We need to add four x and y coordinates to plot a rectangle for each gene now
    
    # Each gene set is on a separate row in the plot. Let's make it so that the
    # most significant result is at the top
    plot_data$y1 <- order(unique(plot_data$gsa_pval), decreasing = TRUE)[match(plot_data$gsa_pval, unique(plot_data$gsa_pval))]
    plot_data <- mutate(plot_data, up_down = case_when(log2FC > 0 ~ "up",
                                            TRUE ~ "down")) %>%
        mutate(y2 = y1 + 1) %>%
        group_by(gene_set, up_down) %>%
        arrange(log2FC) %>%
        mutate(x1=row_number() -1) %>%
        mutate(x1=case_when(up_down=="down" ~ -1*rev(x1),
                            TRUE ~ x1)) %>%
        mutate(x2=case_when(up_down=="down" ~ x1 - 1,
                            TRUE ~ x1 + 1)) %>%
        ungroup()
    
    # Pivot to get one column with x values another with y values
    plot_data <- plot_data %>%
        pivot_longer(cols = c(x1,x2),
                     values_to = 'x',
                     names_to = "x.cord")%>%
        pivot_longer(cols = c(y1,y2),
                     values_to = 'y',
                     names_to = "y.cord") %>%
        select(-x.cord, -y.cord)
    
    # Add group id
    plot_data <- plot_data %>%
        mutate(id = str_glue('{gene_set}-{gene}'))
    
    # Change order of y cords to plot rectangles (not bowties)
    plot_data <- plot_data %>%
        group_by(gene_set, gene) %>%
        mutate(y=y[c(2,1,3,4)])
    
    # Get subset for plotting gene set labels
    text.labels <- plot_data %>%
        group_by(gene_set) %>%
        slice_max(order_by = x, n=1) %>%
        slice_max(order_by = y, n=1) %>%
        mutate(y=y-0.5,
               x=x+1) %>%
        select(gene_set,x,y)
        
       # Create plot
    p <- plot_data %>% 
        ggplot(aes(x = x, y = y)) +
        geom_polygon(aes(fill = log2FC, group = id)) +
        theme_classic() +
        scale_fill_gradient2(low = "blue", mid="grey", high = "red", midpoint = 0) +
        geom_vline(xintercept = 0) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.line.y=element_blank()) +
        scale_x_continuous(labels = abs)
    
    # Add the gene set labels
    p <- p + geom_text(data = text.labels, aes(x=x, y=y, label=gene_set))
    
    # Add text beneath the x axis
    p <- add_sub(p, "# genes down-regulated", x = 0, hjust = 0)
    p <- add_sub(p, "# genes up-regulated", x = 1, hjust = 1, vjust = 0, y = 1)
    
    return(ggdraw(p))
    }

pdf(str_glue("{fig_dir}/{sample_name}_gsa_barplot.pdf"), width=11.6, height=8.2)
gsa_plot()
dev.off()



#Last thing.  Get the top and bottom 10 FC genes for each pathway
#Get the list of pathways in the final table

gene.set.list = gene.sets
fold.change.results = fold.change
pval.results = pvals

    plot_data <- map_df(.x = 1:10, .f = function(i) expand.grid(unlist(gene.set.list[[i]]), names(gene.set.list)[i])) %>%
        rename(gene = Var1, gene_set = Var2) %>%
        merge(fold.change.results, by="gene") %>%
        merge(pval.results, by="gene_set")

table(plot_data$gene_set)

#Subset the table based on the first one

test<-subset(plot_data, gene_set == "Degradation of the extracellular matrix")
#Get top 10 upregulated genes
test<-test[order(-test$log2FC), ]
top<-test[1:10, ]
test<-test[order(test$log2FC), ]
bottom<-test[1:10, ]

#bind together
highFCgenes<-rbind(top, bottom)
highFCgenes<-highFCgenes[order(highFCgenes$log2FC), ]

#Save out
write.csv(highFCgenes, file = str_glue("{gene_data}/{sample_name}_Degradation of the extracellular matrix_highFCgenes.csv"))


#Need to finish combined.
 Regulation of expression of SLITs and ROBOs 
                                                                                                                 76 
                                                                                                        Translation 
                                                                                                                 98 
                                                                                        Signaling by ROBO receptors 
                                                                                                                 81 
Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins. 
                                                                                                                 59 
                                                     The citric acid (TCA) cycle and respiratory electron transport 
                                                                                                                 70 
                                                        SRP-dependent cotranslational protein targeting to membrane 
                                                                                                                 53 
                                                                                     Respiratory electron transport 
                                                                                                                 49 
                                                                 GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2 
                                                                                                                 34 
                                                                        Regulation of ornithine decarboxylase (ODC) 
                                                                                                                 33 
                                                                          Metabolism of amino acids and derivatives 
                                                                                                                 97 