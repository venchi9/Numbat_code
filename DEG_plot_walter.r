# Load R packages ---------------------------------------------------------

library(tidyverse)
library(cowplot)


# Load data ---------------------------------------------------------------

set.seed(123)

# Simulate list of gene sets
gene.sets <- list(geneset1 = letters[1:5],
                  geneset2 = letters[10:26],
                  geneset3 = sample(letters,15),
                  geneset4 = sample(letters,10),
                  geneset5 = sample(letters,15))

# Simulate log fold change values
fold.change <- data.frame(gene = letters,
                          log2FC = runif(length(letters), min = -2, max = 2))

# Get p-values for each gene set
pvals <- data.frame(gene_set = c("geneset1","geneset2","geneset3","geneset4","geneset5"),
                    gsa_pval = c(0.00001, 0.0001, 0.001, 0.01, 0.1))



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
    plot_data <- map_df(.x = 1:5, .f = function(i) expand.grid(unlist(gene.set.list[[i]]), names(gene.set.list)[i])) %>%
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

gsa_plot()