
#Plots are now fixed.~

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

#Samples
sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"

#Comparisons
comparison_name<-"LN_v_primary4"
comparison_name<-"LN_v_primary6"
comparison_name<-"LN2_vP2"
comparison_name<-"LN3_v_P3"
comparison_name<-"LN23_v_P23"
comparison_name<-"LN3_v_primary"
comparison_name<-"LN4_v_primary"
comparison_name<-"LN_v_primary"


#Step one.  Organise data for gene.set.list:  a named list of the genes in each gene set
#Load up the reactome pathways.
geneset<-read.csv(file = str_glue("{data_dir}/{sample_name}_{comparison_name}_reactome_enriched_pathways.csv"))

#Make vectors for the top 10 genelists
x<-geneset[1:10, ]

#Walter's solution
gene.sets <- strsplit(as.character(x$geneID), "/")
names(gene.sets) <- x$Description

#Now get a list of the log2FC for each gene in our data
FC<-read.csv(file = str_glue("{data_dir}/{sample_name}_DEG_analysis_{comparison_name}.csv"))


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



gsa_plot <- function(gene.set.list = gene.sets,
                     fold.change.results = fold.change,
                     pval.results = pvals){
  
  # Fix bug caused by duplicated values 
  pval.results$y1 <- order(pval.results$gsa_pval, decreasing = TRUE)
  
  # Change inputs into a single data frame for plotting
  plot_data <- map_df(.x = 1:nrow(pval.results), .f = function(i) expand.grid(unlist(gene.set.list[[i]]), names(gene.set.list)[i])) %>%
    rename(gene = Var1, gene_set = Var2) %>%
    merge(fold.change.results, by="gene") %>%
    merge(pval.results, by="gene_set")
  
  # We need to add four x and y coordinates to plot a rectangle for each gene now
  # Each gene set is on a separate row in the plot. Let's make it so that the
  # most significant result is at the top
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
    scale_fill_gradientn(colours = c("darkblue","blue","grey", "red", "darkred"), limits =c(-3, 2.5), breaks =c(-2, -1,0,1, 2)) +
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
  p <- cowplot::add_sub(p, "# genes down-regulated", x = 0, hjust = 0)
  p <- cowplot::add_sub(p, "# genes up-regulated", x = 1, hjust = 1, vjust = 0, y = 1)
  
  return(cowplot::ggdraw(p))
}

#Plot
pdf(str_glue("{fig_dir}/{sample_name}_{comparison_name}_gsa_barplot.pdf"), width=11.6, height=8.2)
gsa_plot()
dev.off()