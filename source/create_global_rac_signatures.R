setwd("/data/ruoffcj/projects/drug_treatment/")
library(Seurat)
library(presto)
library(corto)
library(tidyverse)

cell_lines <- c("A549","K562","MCF7")

# RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17)) # OR > 1.5
RACs <- list(c(9,12,13,14,16,18),c(4,5,9,11),c(5,8,13,17)) # OR > 2
names(RACs) <- cell_lines

# Create Global RAC Signatures for each cell line
rac_signatures <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  #read in DR signature scores and set active cells
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  Idents(data) <- data$rac
  
  #Find global RAC markers for current cell line
  de_res <- FindAllMarkers(data)
  
  saveRDS(de_res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_global_rac_de.rds"))
  
  curr_rac_up_genes <- de_res %>% 
    filter(cluster == "rac" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  rac_signatures[[curr_cell_line]] <- curr_rac_up_genes[1:200]
}

saveRDS(rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_signatures.rds")