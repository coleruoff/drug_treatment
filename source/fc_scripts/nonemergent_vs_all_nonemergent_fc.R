library(Seurat)
library(tidyverse)

cell_lines <- c("A549", "K562","MCF7")

for(curr_cell_line in cell_lines){
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
  
  emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
  emergent <- emergent[[curr_cell_line]]
  non_emergent <- levels(data)[-emergent]
  
  clusters_of_interest <- non_emergent
  
  data <- data[,data$Cluster %in% clusters_of_interest]
  
  curr_list <- list()
  
  for(curr_cluster in clusters_of_interest){
    
    curr_fc_results <- FoldChange((data[,data$Cluster %in% non_emergent]), ident.1=curr_cluster)
    
    curr_list <- append(curr_list, list(curr_fc_results))
  }
  
  names(curr_list) <- paste0(curr_cell_line, "_cluster", clusters_of_interest, "_fc_results")
  
  saveRDS(curr_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_nonemergent_vs_all_nonemergent_fc.rds"))
}


