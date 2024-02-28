source("/data/ruoffcj/projects/drug_treatment/source/cole_functions.R")
source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

cell_lines <- c("A549","K562","MCF7")

################################################################################

# Components of RAC supercluster
supercluster_components <- list()
supercluster_components[["A549"]] <- c(4,9,13)
supercluster_components[["K562"]] <- c(5)
supercluster_components[["MCF7"]] <-  c(5,8,17)

all_data <- list()
all_data[["A549"]] <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/A549_processed_filtered.rds")
all_data[["K562"]] <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/K562_processed_filtered.rds")
all_data[["MCF7"]] <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/MCF7_processed_filtered.rds")

################################################################################

component_signatures <- list()

for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  #Read in cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  # 
  # #read in DR signature scores and set active cells
  # scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  # threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  # active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  # clusters_of_interest <- RACs[[curr_cell_line]]
  # 
  # #Add metadata for RAC and Cell Group
  # data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  # data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  # data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  # 
  # Idents(data) <- data$Cluster
  
  # de_results <- FindAllMarkers(data)
  # saveRDS(de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))
  
  de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))
  
  curr_rac_signatures <- list()
  for(i in levels(data)){
    
    curr_signature <- de_results %>%
      filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == i) %>%
      arrange(desc(avg_log2FC)) %>%
      pull(gene)
    
    if(length(curr_signature) > 200){
      curr_signature <- curr_signature[1:200]
    }
    
    curr_rac_signatures <- append(curr_rac_signatures, list(curr_signature))
    
  }
  
  names(curr_rac_signatures) <- paste0(curr_cell_line, "_", levels(data))
  curr_rac_signatures <- curr_rac_signatures[supercluster_components[[curr_cell_line]]]
  
  component_signatures <- append(component_signatures,curr_rac_signatures)
}


all_control_signatures <- list()
for(i in 1:500){
  cat(i, "\n")
  ##############
  # Find supercluster 1 control signature
  
  control_supercluster_signatures <- component_signatures
  for(curr_cell_line in cell_lines){
    
    cat(curr_cell_line, "\n")
    
    curr_signatures <- component_signatures[grepl(curr_cell_line, names(component_signatures))]
    curr_signatures_names <- names(component_signatures)[grepl(curr_cell_line, names(component_signatures))]
    
    
    control_supercluster_signatures[curr_signatures_names] <- find_control_gene_sets(all_data[[curr_cell_line]][["RNA"]]$counts, curr_signatures, num_bins=10)
  }
  
  
  control_rac_supercluster_signatures <- list("supercluster_signature" = find_consensus_geneset(control_supercluster_signatures,2))
  
  
  all_control_signatures <- append(all_control_signatures, (control_rac_supercluster_signatures))
}


saveRDS(all_control_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/control_rac_supercluster_signatures.rds")



