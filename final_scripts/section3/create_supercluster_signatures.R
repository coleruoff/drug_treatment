args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
source("final_scripts/drug_treatment_functions.R")
library(tidyverse)
library(Seurat)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#################################################################################

cell_lines <- c("A549","K562","MCF7")

num_clusters <- c(19,12,18)
names(num_clusters) <- cell_lines

all_signatures <- list()
for(curr_cell_line in cell_lines){

  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  curr_signatures <- list()
  for(curr_cluster in 1:num_clusters[curr_cell_line]){

    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de.rds"))

    curr_cluster_signature <- de_res %>% 
      filter(cluster == curr_cluster & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene) 
    

    curr_signatures[[curr_cluster]] <- curr_cluster_signature#[1:800]

  }

  all_signatures[[curr_cell_line]] <- curr_signatures
}
################################################################################

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

################################################################################

supercluster1_signatures <- c(list(all_signatures[["A549"]][[supercluster_components[[1]][["A549"]]]]),
                              list(all_signatures[["K562"]][[supercluster_components[[1]][["K562"]]]]),
                              list(all_signatures[["MCF7"]][[supercluster_components[[1]][["MCF7"]]]]))

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_signatures,2))

##############

supercluster2_signatures <- c(list(all_signatures[["A549"]][[supercluster_components[[2]][["A549"]]]]),
                              list(all_signatures[["K562"]][[supercluster_components[[2]][["K562"]]]]),
                              list(all_signatures[["MCF7"]][[supercluster_components[[2]][["MCF7"]]]]))

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_signatures,2))

supercluster_signatures <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_signatures, paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

################################################################################
################################################################################
# Downregulated genes

all_signatures <- list()
for(curr_cell_line in cell_lines){
  
  # de_res <- readRDS(paste0(dataDirectory, "de_results/",curr_cell_line,"_cluster_de_MAST.rds"))
  # de_res <- readRDS(paste0(dataDirectory, "de_results/",curr_cell_line,"_cluster_de.rds"))
  
  # curr_clusters <- unique(de_res$cluster)
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  curr_signatures <- list()
  for(curr_cluster in 1:num_clusters[curr_cell_line]){
    
    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de_MAST.rds"))
    
    curr_cluster_signature <- de_res %>%
      filter(cluster == curr_cluster & avg_log2FC < 0 & p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      pull(gene)
    
    curr_signatures[[curr_cluster]] <- curr_cluster_signature#[1:800]
    
  }
  
  all_signatures[[curr_cell_line]] <- curr_signatures
  
}

################################################################################

supercluster1_signatures <- c(list(all_signatures[["A549"]][[9]]),list(all_signatures[["K562"]][[5]]),list(all_signatures[["MCF7"]][[8]]))

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_signatures,3))

##############

supercluster2_signatures <- c(list(all_signatures[["A549"]][[14]]),list(all_signatures[["K562"]][[9]]),list(all_signatures[["MCF7"]][[13]]))

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_signatures,3))

supercluster_signatures <- c(supercluster1_signature,supercluster2_signature)

saveRDS(supercluster_signatures, paste0(dataDirectory, "genesets/rac_supercluster_down_signatures.rds"))

