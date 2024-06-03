library(ggpubr)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


# curr_cell_line <- "A549"
# ident_to_use <- "Cluster"
# geneset_to_use <- "raj_watermelon_resistance_signature"


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

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")

df <- create_df_for_boxplots(data, "Cluster", "pan_stem_genes")

plot_title <- paste0("Resistance Signature AUCell Scores in ", curr_cell_line, " Clusters")

p <- ggboxplot(df, x = "Cluster", y = "curr_geneset", fill="group_median_score")+
  scale_fill_gradient(low="white", high="red", name="Median Cluster Score")+
  ggtitle(plot_title)+
  xlab("Clusters")+
  ylab("Score")+
  theme(legend.position="right",
        title = element_text(size=28, face="bold"),
        axis.text = element_text(size=30),
        legend.text = element_text(size=14),
        legend.title = element_text(size=12))


# png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/resistance_score_boxplots/", curr_cell_line,"_boxplots.png"),
    # width=1500, height=500)
p

# dev.off()