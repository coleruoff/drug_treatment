library(Seurat)
library(tidyverse)
set.seed(42)

cell_lines <- c("A549","K562","MCF7")
geneset_to_use <- "raj_watermelon_resistance_signature"

curr_cell_line <- "K562"

x <- 10
all_active_subpops <- list()

for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  # data <- data[,data$treatment_stage == "post"]
  
  all_clusters <- sort(as.numeric(unique(data$Cluster)))
  num_clusters <- length(all_clusters)
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
  
  active_cells <- rownames(scores)[scores>threshold$threshold]
  
  data <- AddMetaData(data, metadata = scores,col.name = "resistant_score")
  
  data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cells, "active","inactive"),col.name = "active")
  
  subpops <- c()
  differences <- c()
  i <- 4
  for(i in 1:num_clusters){
    cat(i,"\n")

    curr_cluster_names <- colnames(data)[data$Cluster == i]
    
    curr_cluster_active_cells <- curr_cluster_names[curr_cluster_names %in% active_cells]
    
    if(length(curr_cluster_active_cells) > 0){  
      curr_PCs <- data@reductions$PCA@cell.embeddings[curr_cluster_active_cells,1:20]
      
      active_dist <- as.matrix(dist(curr_PCs))
      
      random_p_vals <- c()
      for(j in 1:500){
        cat(j,"\n")
        random_cells <- sample(curr_cluster_names, length(curr_cluster_active_cells))
        
        curr_PCs <- data@reductions$PCA@cell.embeddings[random_cells,1:20]
        
        random_dist <- as.matrix(dist(curr_PCs))
        
        active_dist <- as.vector(active_dist)
        random_dist <- as.vector(random_dist)
        
        wilcox_res <- wilcox.test(active_dist, random_dist)
        
        random_p_vals <- append(random_p_vals, wilcox_res$p.value)
      }
      
      sum(random_p_vals > .05)/length(random_p_vals)
      
      active_dist <- as.vector(active_dist)
      random_dist <- as.vector(random_dist)
      
      wilcox_res <- wilcox.test(active_dist, random_dist)
      
      if(wilcox_res$p.value < 0.05 & median(active_dist) < median(random_dist)){
        subpops <- append(subpops, i)
        differences <- append(differences, (median(active_dist) - median(random_dist)))
        
      }
    }
  }
  all_active_subpops <- append(all_active_subpops, list(subpops))
}



names(all_active_subpops) <- cell_lines
saveRDS(all_active_subpops, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/all_active_subpops.rds")

