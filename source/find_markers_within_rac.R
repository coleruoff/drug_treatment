library(Seurat)
library(tidyverse)
library(org.Hs.eg.db)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')

args <- commandArgs(trailingOnly=TRUE)
curr_cell_line <- args[1]

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")

Idents(data) <- data$cell_group

all_signatures <- list()

for(curr_cluster in clusters_of_interest){
  cat(curr_cluster, "\n")
  
  de_res <- FindAllMarkers(data[,data$Cluster == curr_cluster])
  
  temp <- de_res %>% 
    filter(p_val_adj < 0.05 & cluster == "1" & avg_log2FC > .5) %>% 
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene,avg_log2FC)
  
  # curr_signature <- temp$avg_log2FC
  # names(curr_signature) <- temp$gene
  
  curr_signature <- temp$gene
  
  all_signatures <- append(all_signatures, list(curr_signature))
}

names(all_signatures) <- clusters_of_interest

saveRDS(all_signatures, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_rac_within_cluster_de_signatures.rds"))
