setwd("/data/ruoffcj/projects/drug_treatment/")
library(ggpubr)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

create_df_for_boxplots <- function(data, ident_to_use, geneset_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(geneset_to_use, "\n")
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
  
  # Create list of lists containing active cells for each geneset
  all_active_cells <- list()
  
  curr_geneset <- colnames(scores)[1]
  # for(curr_geneset in colnames(scores)){
  
  cat(curr_geneset, "\n")
  active_cells <- rownames(scores)[scores[,curr_geneset] > threshold$threshold[threshold$gene_set == curr_geneset]]
  
  # all_active_cells <- append(all_active_cells, list(curr_active_cells))
  
  data <- AddMetaData(data, metadata = scores[,curr_geneset], col.name = "curr_geneset")
  
  
  df <- data@meta.data %>% 
    dplyr::select(ident_to_use, curr_geneset, dose, treatment_stage) %>% 
    group_by(eval(parse(text=ident_to_use))) %>% 
    mutate("group_median_score" = median(curr_geneset))
  # }
  
  return(df)    
}

###################################################################################

curr_cell_line <- "MCF7"

cell_lines <- c("A549","K562","MCF7")

plots <- list()

for(curr_cell_line in cell_lines){
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
  
  df <- create_df_for_boxplots(data, "Cluster", "raj_watermelon_resistance_signature")
  
  # plot_title <- paste0("Resistance Signature AUCell Scores in ", curr_cell_line, " Clusters")
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
  
  plots <- append(plots,list(p))
}


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_1b.png"),
     width=20, height=20, units = "in", res = 300)

figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = textGrob("Score", rot = 90, vjust = 1, gp = gpar(fontsize=35, fontface="bold")),
                bottom = textGrob("Clusters", gp = gpar(fontsize=35, fontface="bold")),
                top=textGrob("Resistance Signature AUCell Scores in Cell Line Clusters", gp = gpar(fontsize=40, fontface="bold")))

print(p)

dev.off()
