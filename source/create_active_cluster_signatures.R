library(tidyverse)
source("source/cole_functions.R")


cell_lines <- c("A549","K562","MCF7")
geneset_to_use <- "raj_resistant_2017"

curr_cell_line <- cell_lines[3]

active_cluster_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_", geneset_to_use, "_all_active_cluster_markers.rds"))


all_clusters <- sort(as.character(unique(active_cluster_de$cluster)))
i <- all_clusters[1]

active_cluster_signatures <- list()

#Create signatures
for(i in all_clusters){
  
  curr_signature <- active_cluster_de %>% 
    filter(cluster == i & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  # names(curr_signature) <- active_cluster_de %>% 
  #   filter(cluster == i & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  #   arrange(desc(avg_log2FC)) %>% 
  #   pull(gene)
  
  
  active_cluster_signatures <- append(active_cluster_signatures, list(curr_signature))
  
}

names(active_cluster_signatures) <- paste0("cluster", all_clusters, "_signature")

active_clusters_names <- names(active_cluster_signatures)[grepl("[0-9]a",names(active_cluster_signatures))]
#Remove shared genes in clusters and the active subpop in that cluster
for(i in active_clusters_names){
  
  
  inactive_component_name <- sub("a", "",i)
  
  shared_genes <- intersect(active_cluster_signatures[[i]], active_cluster_signatures[[inactive_component_name]])
  
  active_cluster_signatures[[i]] <- active_cluster_signatures[[i]][!active_cluster_signatures[[i]] %in% shared_genes]
  
  active_cluster_signatures[[inactive_component_name]] <- active_cluster_signatures[[inactive_component_name]][!active_cluster_signatures[[inactive_component_name]] %in% shared_genes]  
}


jaccard_mat <- calc_jaccard_matrix(active_cluster_signatures, active_cluster_signatures)


Heatmap(jaccard_mat)














