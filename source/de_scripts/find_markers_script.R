library(Seurat)
library(tidyverse)

#Runs FindMarkers on all three sci-Plex cell lines with Trapnell Clusters as Idents

cell_lines <- c("A549","K562","MCF7")

for(curr_cell_line in cell_lines){
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  de_result <- FindAllMarkers(data, logfc.threshold = .25, test.use = "MAST")
  
  saveRDS(de_result, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_cluster_all_markers_de_25.rds"))
}