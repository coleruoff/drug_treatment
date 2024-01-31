library(Seurat)
library(tidyverse)

cell_lines <- c("A549","K562","MCF7")

curr_cell_line <- cell_lines[1]
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))
  
  Idents(data) <- data$treatment_stage
  
  all_clusters <- sort(as.numeric(unique(data$Cluster)))
  
  cluster_ranked_markers <- list()
  i <- 2
  for(i in all_clusters){
    curr_cluster_markers <- FindAllMarkers(data[,data$Cluster == i])
    
    curr_cluster_post_treatment_up <- curr_cluster_markers %>% 
      filter(cluster == "post" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(avg_log2FC)
      
    names(curr_cluster_post_treatment_up) <- curr_cluster_markers %>% 
      filter(cluster == "post" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    curr_cluster_post_treatment_down <- curr_cluster_markers %>% 
      filter(cluster == "post" & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
      arrange(avg_log2FC) %>% 
      pull(avg_log2FC)
    
    names(curr_cluster_post_treatment_down) <- curr_cluster_markers %>% 
      filter(cluster == "post" & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
      arrange(avg_log2FC) %>% 
      pull(gene)
    
    curr_cluster_up_and_down_list <- list(curr_cluster_post_treatment_up,curr_cluster_post_treatment_down)
    
    names(curr_cluster_up_and_down_list) <- c("up","down")
    
    cluster_ranked_markers <- append(cluster_ranked_markers, list(curr_cluster_up_and_down_list))
    
  }
  
  
  saveRDS(post_treatment_cluster_ranked_markers, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_pre_to_post_cluster_markers.rds"))
}

