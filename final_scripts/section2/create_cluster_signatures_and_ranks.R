args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(Seurat)
library(tidyverse)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

num_clusters <- c(19,12,18)
names(num_clusters) <- cell_lines

length_to_use <- 200

all_signatures <- list()
all_ranks <- list()
for(curr_cell_line in cell_lines){
  
  # curr_cell_line <- cell_lines[1]
  
  # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de_MAST.rds"))
  # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  # curr_cell_line_signatures <- list()
  # curr_cell_line_ranks <- list()
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  for(curr_cluster in as.character(1:num_clusters[curr_cell_line])){
    
    # curr_cluster <- 1
    
    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de.rds"))
    
    # Create current cluster gene signature
    curr_signature <- de_res %>% 
      filter(cluster == curr_cluster & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene) 
  
    
    # Trim signature to set length
    if(length(curr_signature) > 200){
      curr_signature <- curr_signature[1:length_to_use] 
    }
    
    all_signatures[[curr_cell_line]][[curr_cluster]] <- curr_signature
    
    # Create current cluster gene ranks
    curr_ranks <- de_res %>% 
      arrange(desc(avg_log2FC)) %>% 
      rownames_to_column() %>% 
      dplyr::select(rowname,avg_log2FC) %>% 
      deframe()
    
    all_ranks[[curr_cell_line]][[curr_cluster]] <- curr_ranks
    
  }
}

saveRDS(all_signatures, paste0(dataDirectory, "genesets/cluster_signatures.rds"))
saveRDS(all_ranks, paste0(dataDirectory, "genesets/cluster_ranks.rds"))

################################################################################
# Down regulated signatures

all_signatures <- list()
all_ranks <- list()
for(curr_cell_line in cell_lines){
  
  # curr_cell_line <- cell_lines[1]
  
  # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de_MAST.rds"))
  # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  # curr_cell_line_signatures <- list()
  # curr_cell_line_ranks <- list()
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  for(curr_cluster in as.character(1:num_clusters[curr_cell_line])){
    
    # curr_cluster <- 1
    
    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de.rds"))
    
    # Create current cluster gene signature
    curr_signature <- de_res %>% 
      filter(cluster == curr_cluster & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
      arrange(avg_log2FC) %>% 
      pull(gene)
    
    # Trim signature to set length
    if(length(curr_signature) > 200){
      curr_signature <- curr_signature[1:length_to_use] 
    }
    all_signatures[[curr_cell_line]][[curr_cluster]] <- curr_signature
    
  }
}

saveRDS(all_signatures, paste0(dataDirectory, "genesets/cluster_down_signatures.rds"))




