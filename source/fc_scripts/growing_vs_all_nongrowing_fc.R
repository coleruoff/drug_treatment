library(Seurat)
library(tidyverse)

cell_lines <- c("A549", "K562","MCF7")

curr_cell_line <- cell_lines[1]

for(curr_cell_line in cell_lines){
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
  
  growing <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/growing_clusters.rds")
  growing <- sort(unique(unlist(growing[[curr_cell_line]])))
  non_growing <- levels(data)[-growing]
  
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% growing, "growing", "non-growing"), col.name = "growing")
  
  Idents(data) <- data$growing
  
  curr_list <- list()
  
  for(curr_growing in growing){
    
    # curr_data <- data[,data$Cluster %in% non_emergent | data$Cluster == curr_emergent]
    
    curr_fc_results <- FoldChange((data[,data$Cluster %in% non_growing | data$Cluster == curr_growing]), ident.1="growing")
    
    curr_list <- append(curr_list, list(curr_fc_results))
    
  }
  
  names(curr_list) <- paste0(curr_cell_line, "_cluster", growing, "_fc_results")
  
  saveRDS(curr_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_growing_vs_all_nongrowing_fc.rds"))
}