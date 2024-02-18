library(Seurat)
library(tidyverse)
library(presto)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")



pre_racs <- c(4,9,12,13)
emergent_racs <- c(14,16,18,19)


curr_cell_line <- cell_lines[1]
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$Cluster %in% pre_racs, "pre_rac", ifelse(data$rac == "rac" & data$Cluster %in% emergent_racs, "emergent_rac", "non-rac")), col.name = "rac_type")
  

  Idents(data) <- data$rac_type
  
  de_results <- FindMarkers(data, ident.1 = "pre_rac",ident.2 = "emergent_rac")
  
  
}




pre_rac_genes <- de_results %>%
  filter(p_val_adj < 0.05, avg_log2FC < 0) %>% 
  arrange(avg_log2FC) %>% 
  rownames()


pre_rac_genes <- list("pre_rac_genes" = pre_rac_genes)


all_dotplots <- genesets_characterization(pre_rac_genes, universe_to_use = cell_line_universes[[curr_cell_line]])




all_dotplots



  
  
  
  
  
  
  
  




