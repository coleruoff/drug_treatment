args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

library(ggpubr)
library(Seurat)
library(tidyverse)
library(grid)

# dataDirectory <- "//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

###################################################################################

cell_lines <- c("A549","K562","MCF7")

plots <- list()

for(curr_cell_line in cell_lines){
  
  # Read in current cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  # Read in resistance signature AUCell scores
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  
  # Add AUCell scores to metadata then subset Cluster and scores for boxplot plotting
  data <- AddMetaData(data, metadata = scores, col.name = "curr_geneset")
  df <- data@meta.data %>% 
    dplyr::select(Cluster, curr_geneset) %>% 
    group_by(Cluster) %>% 
    mutate("group_median_score" = median(curr_geneset))
  
  # Plot current cell line's boxplots
  plot_title <- curr_cell_line
  p <- ggboxplot(df, x = "Cluster", y = "curr_geneset", fill="group_median_score")+
    scale_fill_gradient(low="white", high="orange", name="      Median\nCluster Score")+
    ggtitle(plot_title)+
    xlab("")+
    ylab("")+
    theme(legend.position="right",
          title = element_text(size=24, face="bold"),
          axis.text = element_text(size=30),
          legend.text = element_text(size=24),
          legend.title = element_text(size=26),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1,"cm"))
  
  # Add current cell line boxplots to final list of plots
  plots <- append(plots,list(p))
}

# Arrange each cell line plot, annotate figure, and plot
figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = textGrob("Score", rot = 90, vjust = 1, gp = gpar(fontsize=35, fontface="bold")),
                bottom = textGrob("Clusters", gp = gpar(fontsize=35, fontface="bold")),
                top=textGrob("Resistance Signature AUCell Scores in Cell Line Clusters", gp = gpar(fontsize=40, fontface="bold")))

png(paste0(plotDirectory,"figure_1b.png"),
    width=20, height=20, units = "in", res = 300)

print(p)

dev.off()
