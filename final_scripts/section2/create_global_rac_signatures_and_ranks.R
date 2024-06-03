args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(Seurat)
library(tidyverse)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

################################################################################
cell_lines <- c("A549","K562","MCF7")

length_to_use <- 200

all_signatures <- list()
all_ranks <- list()
for(curr_cell_line in cell_lines){
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_global_rac_de_MAST.rds"))
  # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_global_rac_de.rds"))
  
  curr_cell_line_signatures <- list()
  curr_cell_line_ranks <- list()
  
  # Create global RAC gene signature
  curr_signature <- de_res %>% 
    filter(cluster == "rac" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  # Trim signature to set length
  curr_signature <- curr_signature[1:length_to_use] 
  
  all_signatures[[curr_cell_line]] <- curr_signature
  
  # Create global RAC gene ranks
  curr_ranks <- de_res %>% 
    filter(cluster == "rac") %>% 
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene,avg_log2FC) %>% 
    deframe()
  
  # Trim ranks to set length
  all_ranks[[curr_cell_line]] <- curr_ranks
  
}

saveRDS(all_signatures, paste0(dataDirectory, "genesets/global_rac_signatures.rds"))
saveRDS(all_ranks, paste0(dataDirectory, "genesets/global_rac_ranks.rds"))




