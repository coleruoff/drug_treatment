library(tidyverse)


dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

cell_lines <- c("A549","K562","MCF7")
curr_cell_line <- cell_lines[1]

trapnell_de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))

all_clusters <- as.numeric(sort(unique(trapnell_de_res$cluster)))

signature_length <- 500

cluster_signatures <- list()

for(i in all_clusters){
  
  curr_markers <- trapnell_de_res %>% 
    filter(cluster == i & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(avg_log2FC) %>% 
    pull(gene)
  
  if(length(curr_markers) > signature_length){
    curr_markers <- curr_markers[1:signature_length]  
  }
  
  cluster_signatures <- append(cluster_signatures, list(curr_markers))

}


names(cluster_signatures) <- paste0(curr_cell_line, "_cluster", all_clusters,"_signature")

saveRDS(cluster_signatures, paste0(dataDirectory, "genesets/", curr_cell_line, "_cluster_signatures_", signature_length, ".rds"))


