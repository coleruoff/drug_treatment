library(Seurat)
library(tidyverse)

cell_lines <- c("A549", "K562","MCF7")

for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
  
  growing <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/growing_clusters.rds")
  growing <- sort(unique(unlist(growing[[curr_cell_line]])))
  non_growing <- levels(data)[-growing]
  
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% growing, "growing", "non-growing"), col.name = "growing")
  
  Idents(data) <- data$growing
  
  curr_list <- list()
  
  for(curr_growing in growing){
    
    curr_de_results <- FindMarkers((data[,data$Cluster %in% non_growing | data$Cluster == curr_growing]), ident.1="growing")
    
    curr_list <- append(curr_list, list(curr_de_results))
    
  }
  
  names(curr_list) <- paste0(curr_cell_line, "_cluster", growing, "_de_results")
  
  saveRDS(curr_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_growing_vs_nongrowing_de.rds"))
}