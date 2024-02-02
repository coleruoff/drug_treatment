library(Seurat)
library(tidyverse)
library(presto)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

#################################################################################

# rac_signatures <- list()
rac_type_signatures <- list()
rac_type1_signatures <- list()
rac_type2_signatures <- list()
intra_rac_signatures <- list()

rac_type_down_signatures <- list()
rac_type1_down_signatures <- list()
rac_type2_down_signatures <- list()
intra_rac_down_signatures <- list()


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
  
  #################################################################################
  # Create Gene Signatures
  #################################################################################
  
  #################################################################################
  # RAC Signatures

  # Idents(data) <- data$Cluster
  #
  # # de_results <- FindAllMarkers(data)
  # # saveRDS(de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))
  #
  # de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))
  #
  # curr_rac_signatures <- list()
  # for(i in levels(data)){
  #
  #   curr_signature <- de_results %>%
  #     filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == i) %>%
  #     arrange(desc(avg_log2FC)) %>%
  #     pull(gene)
  #
  #   if(length(curr_signature) > 200){
  #     curr_signature <- curr_signature[1:200]
  #   }
  #
  #   curr_rac_signatures <- append(curr_rac_signatures, list(curr_signature))
  #
  # }
  #
  # names(curr_rac_signatures) <- paste0(curr_cell_line, "_cluster",levels(data),"_signature")
  # curr_rac_signatures <- curr_rac_signatures[clusters_of_interest]

  #################################################################################
  # Global RAC Type 1 Signature and Global RAC Type 2 Signature
  Idents(data) <- data$cell_group

  # de_results <- FindAllMarkers(data)
  # saveRDS(de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))

  de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))

  curr_rac_type_signatures <- list()
  curr_rac_type_down_signatures <- list()
  
  for(i in levels(data)){

    curr_signature <- de_results %>%
      filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == i) %>%
      arrange(desc(avg_log2FC)) %>%
      pull(gene)
    
    curr_down_signature <- de_results %>%
      filter(p_val_adj < 0.05 & avg_log2FC < 0 & cluster == i) %>%
      arrange(avg_log2FC) %>%
      pull(gene)
    

    if(length(curr_signature) > 200){
      curr_signature <- curr_signature[1:200]
    }
    
    if(length(curr_down_signature) > 200){
      curr_down_signature <- curr_down_signature[1:200]
    }

    
    curr_rac_type_signatures <- append(curr_rac_type_signatures, list(curr_signature))
    curr_rac_type_down_signatures <- append(curr_rac_type_down_signatures, list(curr_down_signature))

  }

  names(curr_rac_type_signatures) <- paste0(curr_cell_line, "_global_rac_type", levels(data), "_signature")
  names(curr_rac_type_down_signatures) <- paste0(curr_cell_line, "_global_rac_type", levels(data), "_down_signature")

  curr_rac_type_signatures <- curr_rac_type_signatures[!grepl("(0)",names(curr_rac_type_signatures))]
  curr_rac_type_down_signatures <- curr_rac_type_down_signatures[!grepl("(0)",names(curr_rac_type_down_signatures))]
  
  #################################################################################
  # Individual RAC Type 1 Signature and Individual RAC Type 2 Signature
  Idents(data) <- data$cell_cluster_group
  
  # de_results <- FindAllMarkers(data)
  # saveRDS(de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_cluster_group_de.rds"))
  
  de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_cluster_group_de.rds"))
  
  curr_rac_subcluster_signatures <- list()
  curr_rac_subcluster_down_signatures <- list()
  for(i in levels(data)){
    
    curr_signature <- de_results %>% 
      filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == i) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    curr_down_signature <- de_results %>% 
      filter(p_val_adj < 0.05 & avg_log2FC < 0 & cluster == i) %>% 
      arrange(avg_log2FC) %>% 
      pull(gene)
    
    if(length(curr_signature) > 200){
      curr_signature <- curr_signature[1:200]
    }
    
    if(length(curr_down_signature) > 200){
      curr_down_signature <- curr_down_signature[1:200]
    }
    
    curr_rac_subcluster_signatures <- append(curr_rac_subcluster_signatures, list(curr_signature))
    curr_rac_subcluster_down_signatures <- append(curr_rac_subcluster_down_signatures, list(curr_down_signature))
    
  }
  
  names(curr_rac_subcluster_signatures) <- paste0(curr_cell_line, "_rac_subcluster", levels(data), "_signature")
  names(curr_rac_subcluster_down_signatures) <- paste0(curr_cell_line, "_rac_subcluster", levels(data), "_down_signature")
  
  curr_rac_type1_signatures <- curr_rac_subcluster_signatures[grepl("(_1)",names(curr_rac_subcluster_signatures))]
  curr_rac_type1_down_signatures <- curr_rac_subcluster_down_signatures[grepl("(_1)",names(curr_rac_subcluster_down_signatures))]
  
  curr_rac_type2_signatures <- curr_rac_subcluster_signatures[grepl("(_2)",names(curr_rac_subcluster_signatures))]
  curr_rac_type2_down_signatures <- curr_rac_subcluster_down_signatures[grepl("(_2)",names(curr_rac_subcluster_down_signatures))]
  
  
  #################################################################################
  # Global RAC Type 1 vs RAC Type 2 Signature
  
  Idents(data) <- data$cell_group
  
  temp_data <- data[,data$cell_group != 0]
  
  de_results <- FindAllMarkers(temp_data)
  
  saveRDS(de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_global_intra_rac_de.rds"))
  
  de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_global_intra_rac_de.rds"))
  
  curr_signature <- de_results %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == 1) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  curr_down_signature <- de_results %>% 
    filter(p_val_adj < 0.05 & avg_log2FC < 0 & cluster == 1) %>% 
    arrange(avg_log2FC) %>% 
    pull(gene)
  
  if(length(curr_signature) > 200){
    curr_signature <- curr_signature[1:200]
  }
  
  if(length(curr_down_signature) > 200){
    curr_down_signature <- curr_down_signature[1:200]
  }
  
  curr_intra_rac_signature <- list(curr_signature)
  curr_intra_rac_down_signature <- list(curr_down_signature)
  
  names(curr_intra_rac_signature) <- paste0(curr_cell_line,"_global_intra_rac_signature")
  names(curr_intra_rac_down_signature) <- paste0(curr_cell_line,"_global_intra_rac_down_signature")
  

  #################################################################################
  # rac_signatures <- append(rac_signatures,curr_rac_signatures)
  rac_type_signatures <- append(rac_type_signatures,curr_rac_type_signatures)
  rac_type1_signatures <- append(rac_type1_signatures,curr_rac_type1_signatures)
  rac_type2_signatures <- append(rac_type2_signatures,curr_rac_type2_signatures)
  intra_rac_signatures <- append(intra_rac_signatures,curr_intra_rac_signature)
  
  rac_type_down_signatures <- append(rac_type_down_signatures,curr_rac_type_down_signatures)
  rac_type1_down_signatures <- append(rac_type1_down_signatures,curr_rac_type1_down_signatures)
  rac_type2_down_signatures <- append(rac_type2_down_signatures,curr_rac_type2_down_signatures)
  intra_rac_down_signatures <- append(intra_rac_down_signatures,curr_intra_rac_down_signature)
}



# saveRDS(rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")
saveRDS(rac_type_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
saveRDS(rac_type1_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")
saveRDS(rac_type2_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type2_signatures.rds")
saveRDS(intra_rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_intra_rac_signatures")

saveRDS(rac_type_down_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_down_signatures.rds")
saveRDS(rac_type1_down_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_down_signatures.rds")
saveRDS(rac_type2_down_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type2_down_signatures.rds")
saveRDS(intra_rac_down_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_intra_rac_down_signatures")


