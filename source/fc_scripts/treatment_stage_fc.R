library(Seurat)
library(tidyverse)

cell_lines <- c("A549", "K562","MCF7")

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
  emergent <- emergent[[curr_cell_line]]
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
  
  curr_list <- list()
  
  all_clusters <- as.numeric(sort(unique(data$Cluster)))
  
  for(curr_cluster in all_clusters){
    
    if(curr_cluster %in% emergent){
      
      data <- data[,data$Cluster == curr_cluster & data$treatment_stage == "pre"]
      
      Idents(data) <- data$treatment_stage
      
      curr_fc_results <- FoldChange(data, ident.1="post")
      
      
    } else {
      data <- data[,data$Cluster == curr_cluster]
      
      Idents(data) <- data$treatment_stage
      
      curr_fc_results <- FoldChange(data, ident.1="post")
    }
    
    
    curr_list <- append(curr_list, list(curr_fc_results))
    
  }
  
  names(curr_list) <- paste0(curr_cell_line, "_cluster", all_clusters, "_fc_results")
  
  saveRDS(curr_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_all_clusters_treatment_stage_fc.rds"))
}
