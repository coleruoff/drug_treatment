
cell_lines <- c("A549","K562","MCF7")

curr_cell_line <- cell_lines[1]

pre_treatment_cluster_ranked_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_pre_treatment_cluster_markers.rds"))

for(i in 1:length(pre_treatment_cluster_ranked_markers)){
  
  pre_treatment_cluster_ranked_markers[[i]] <- names(pre_treatment_cluster_ranked_markers[[i]])
  
}


post_treatment_cluster_ranked_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_post_treatment_cluster_markers.rds"))

for(i in 1:length(post_treatment_cluster_ranked_markers)){
  
  post_treatment_cluster_ranked_markers[[i]] <- names(post_treatment_cluster_ranked_markers[[i]])
  
}

jaccard_mat <- calc_jaccard_matrix(pre_treatment_cluster_ranked_markers, post_treatment_cluster_ranked_markers)

Heatmap(jaccard_mat)
