library(Seurat)
library(tidyverse)

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]


data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", paste0(data$Cluster, "_resistant"),data$Cluster), col.name = "resistant_cluster")

intracluster_active_de_signatures <- list()
for(curr_cluster in clusters_of_interest){
  
  curr_data <- data[,data$Cluster == curr_cluster]
  
  Idents(curr_data) <- curr_data$resistant_active
  
  de_res <- FindAllMarkers(curr_data)
  
  curr_signature <- de_res %>% 
    filter(cluster == "active" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  intracluster_active_de_signatures <- append(intracluster_active_de_signatures, list(curr_signature))
}

names(intracluster_active_de_signatures) <- paste0(clusters_of_interest,"_active_signature")

saveRDS(intracluster_active_de_signatures, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line, "intracluster_active_de_signatures.rds"))






