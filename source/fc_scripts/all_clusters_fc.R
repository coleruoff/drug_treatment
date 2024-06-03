library(Seurat)
library(tidyverse)

cell_lines <- c("A549", "K562","MCF7")

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
  
  # data <- data[,data$pathway == "MAPK"]
  data <- data[,data$treatment_stage == "post"]
  
  curr_list <- list()
  
  all_clusters <- as.numeric(sort(unique(data$Cluster)))
  
  for(curr_cluster in all_clusters){
    
    curr_fc_results <- FoldChange(data, ident.1=curr_cluster)
    
    curr_list <- append(curr_list, list(curr_fc_results))
    
  }
  
  names(curr_list) <- paste0(curr_cell_line, "_cluster", all_clusters, "_fc_results")
  
  saveRDS(curr_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_all_clusters_post_fc.rds"))
}
