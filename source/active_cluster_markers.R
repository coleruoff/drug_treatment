library(Seurat)


cell_lines <- c("A549","K562","MCF7")
all_active_subpops <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/all_active_subpops.rds")
geneset_to_use <- "raj_resistant_2017"
#################################################################################
# Read in Trapnell data
#################################################################################
curr_cell_line <- cell_lines[1]

for(curr_cell_line in cell_lines){
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
  
  active_cells <- rownames(scores)[scores[,geneset_to_use] > threshold]
  
  # curr_active_subpops <- all_active_subpops[[curr_cell_line]]
  # active_metadata <- ifelse(colnames(data) %in% active_cells & data$Cluster %in% curr_active_subpops, paste0(data$Cluster,"a"),data$Cluster)
  
  active_metadata <- ifelse(colnames(data) %in% active_cells, paste0(data$Cluster,"a"),data$Cluster)
  
  data <- AddMetaData(data, metadata = active_metadata, col.name = "active_cluster")
  
  Idents(data) <- data$active_cluster
  
  markers.de <- FindAllMarkers(data)
  
  saveRDS(markers.de, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_", geneset_to_use, "_all_active_cluster_markers.rds"))
}

