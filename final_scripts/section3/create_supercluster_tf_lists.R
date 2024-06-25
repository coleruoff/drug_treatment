library(readxl)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

i <- 2
for(i in 1:2){
  
  curr_tf_activity <- list()
  for(curr_cell_line in cell_lines){
    
    # Get TF activity for all current TFs in each component cell line for current supercluster
    tf_heatmap <- readRDS(paste0(dataDirectory, "processed_data/", curr_cell_line, "_tf_heatmap.rds"))
    
    curr_tf_activity <- append(curr_tf_activity, list(tf_heatmap[,as.character(supercluster_components[[i]][[curr_cell_line]])]))
    
  }
  
  common_tfs <- intersect(intersect(names(curr_tf_activity[[1]]),names(curr_tf_activity[[2]])),names(curr_tf_activity[[3]]))
  
  # Filter to retain common TFs
  for(j in 1:3){
    curr_tf_activity[[j]] <- curr_tf_activity[[j]][names(curr_tf_activity[[j]]) %in% common_tfs]
  }
  
  temp <- data.frame(curr_tf_activity)
  
  # Calculate mean TF activity across cell lines
  mean_sc_tf_activity <- rowMeans(temp)
  
  # mean_sc_tf_activity <- apply(temp, MARGIN = 1, FUN = function(x) min(x))
  
  if(i == 1){
    
    sc1_top_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity > quantile(mean_sc_tf_activity,probs = .75)],decreasing = T)
    
    # sc1_top_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity > 1],decreasing = T)
    
    sc1_top_tfs <-  list("supercluster1" = names(sc1_top_tfs_with_scores))
    
    sc1_bottom_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity < quantile(mean_sc_tf_activity,probs = .25)])
    
    # sc1_bottom_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity < 1])
    
    sc1_bottom_tfs <-  list("supercluster1" = names(sc1_bottom_tfs_with_scores))
    
  } else {
    sc2_top_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity > quantile(mean_sc_tf_activity,probs = .75)],decreasing = T)
    
    # sc2_top_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity > 1],decreasing = T)
    
    sc2_top_tfs <-  list("supercluster1" = names(sc2_top_tfs_with_scores))
    
    sc2_bottom_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity < quantile(mean_sc_tf_activity,probs = .25)])
    # sc2_bottom_tfs_with_scores <- sort(mean_sc_tf_activity[mean_sc_tf_activity < 1])
    
    sc2_bottom_tfs <-  list("supercluster1" = names(sc2_bottom_tfs_with_scores))
    
  }
}


supercluster_tf_list <- c(sc1_top_tfs,sc2_top_tfs)

saveRDS(supercluster_tf_list, paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))

supercluster_tf_list <- c(sc1_bottom_tfs,sc2_bottom_tfs)

saveRDS(supercluster_tf_list, paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))


sc1_top_tfs_with_scores[1:10]

sc2_top_tfs_with_scores[1:10]

