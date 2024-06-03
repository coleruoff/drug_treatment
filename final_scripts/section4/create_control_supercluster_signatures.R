source("source/final_scripts/drug_treatment_functions.R")
source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")
library(tidyverse)
library(Seurat)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"
args = commandArgs(trailingOnly=TRUE)
dataDirectory <- args[1]
plotDirectory <- args[2]

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

#################################################################################

cell_lines <- c("A549","K562","MCF7")

num_clusters <- c(19,12,18)
names(num_clusters) <- cell_lines

all_signatures <- list()
for(curr_cell_line in cell_lines){
  
  curr_signatures <- list()
  for(curr_cluster in 1:num_clusters[curr_cell_line]){
    
    de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de_MAST.rds"))
    
    curr_cluster_signature <- de_res %>%
      filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      rownames()
    
    curr_signatures[[curr_cluster]] <- curr_cluster_signature#[1:800]
    
  }
  
  all_signatures[[curr_cell_line]] <- curr_signatures
}

################################################################################
supercluster1_signatures <- c(list(all_signatures[["A549"]][[9]]),list(all_signatures[["K562"]][[5]]),list(all_signatures[["MCF7"]][[8]]))
names(supercluster1_signatures) <- cell_lines

supercluster2_signatures <- c(list(all_signatures[["A549"]][[14]]),list(all_signatures[["K562"]][[9]]),list(all_signatures[["MCF7"]][[13]]))
names(supercluster2_signatures) <- cell_lines

################################################################################
control_supercluster1_signatures <- list()
control_supercluster2_signatures <- list()

curr_cell_line <- cell_lines[1]

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  curr_signatures <- c(list(supercluster1_signatures[[curr_cell_line]]),list(supercluster2_signatures[[curr_cell_line]]))
  names(curr_signatures) <- c("sc1", "sc2")
  
  curr_cell_line_control_signatures <- find_control_gene_sets(all_data[[curr_cell_line]][["RNA"]]$counts, curr_signatures, num_bins=10)
  
  control_supercluster1_signatures[[curr_cell_line]] <- curr_cell_line_control_signatures$sc1_Control
  control_supercluster2_signatures[[curr_cell_line]] <- curr_cell_line_control_signatures$sc2_Control
}


control_supercluster1_signature <- list("control_supercluster1_signature" = find_consensus_geneset(control_supercluster1_signatures,3))
control_supercluster2_signature <- list("control_supercluster2_signature" = find_consensus_geneset(control_supercluster2_signatures,3))

control_supercluster_signatures <- c(control_supercluster1_signature,control_supercluster2_signature)

saveRDS(control_supercluster_signatures, paste0(dataDirectory, "genesets/rac_supercluster_control_signatures.rds"))




