library(Seurat)
library(tidyverse)

cell_lines <- c("A549","K562","MCF7")


curr_cell_line <- cell_lines[1]
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))
  
  post_treatment_markers <- FindAllMarkers(data[,data$treatment_stage == "post"])
  
  post_treatment_cluster_ranked_markers <- list()
  
  all_clusters <- unique(post_treatment_markers$cluster)
  for(i in all_clusters){
    
    curr_cluster_ranked_markers <- post_treatment_markers %>% 
      filter(cluster == i & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(avg_log2FC)
    
    names(curr_cluster_ranked_markers) <- post_treatment_markers %>% 
      filter(cluster == i & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    
    post_treatment_cluster_ranked_markers <- append(post_treatment_cluster_ranked_markers, list(curr_cluster_ranked_markers))
    
  }
  
  names(post_treatment_cluster_ranked_markers) <- paste0("cluster_",all_clusters)
  
  saveRDS(post_treatment_cluster_ranked_markers, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_post_treatment_cluster_markers.rds"))
}

