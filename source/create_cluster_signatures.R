
cell_lines <- c("A549","K562","MCF7")

################################################################################
# RAC Signatures
################################################################################

all_signatures <- list()

for(curr_cell_line in cell_lines){
  
  curr_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_200.rds"))
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  curr_signatures <- curr_signatures[clusters_of_interest]
  
  for(i in 1:length(curr_signatures)){
    if(length(curr_signatures[[i]]) > 200){
      curr_signatures[[i]] <- curr_signatures[[i]][1:200]
    }
  }
  
  
  all_signatures <- append(all_signatures, curr_signatures)
}

rac_signatures <- all_signatures

saveRDS(rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

################################################################################
# Global active vs inactive signatures
################################################################################

all_signatures <- list()

for(curr_cell_line in cell_lines){
  
  de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_all_resistant_cells_de.rds"))
  
  
  curr_resistant_genes <- de_result %>% 
    filter(cluster == "resistant" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  
  if(length(curr_resistant_genes) > 100){
    curr_resistant_genes <- curr_resistant_genes[1:100]  
  }
  
  
  
  all_signatures <- append(all_signatures, list(curr_resistant_genes))
}

names(all_signatures) <- cell_lines

resistant_signatures <- all_signatures

saveRDS(resistant_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_resistant_signatures.rds")
################################################################################
# Active within RAC signatures
################################################################################

all_signatures <- list()

for(curr_cell_line in cell_lines){
  
  curr_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_rac_within_cluster_de_signatures.rds"))
  
  for(i in 1:length(curr_signatures)){
    if(length(curr_signatures[[i]]) > 100){
      curr_signatures[[i]] <- curr_signatures[[i]][1:100]
    }
  }
  
  names(curr_signatures) <- paste0(curr_cell_line, "_", names(curr_signatures)) 
  
  all_signatures <- append(all_signatures, curr_signatures)
}


within_rac_signatures <- all_signatures

saveRDS(within_rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/within_rac_signatures.rds")

################################################################################
# Curated active subpopulaton signatures
################################################################################

all_active_signatures <- list()
for(curr_cell_line in cell_lines){
  de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_resistant_clusters_de.rds"))
  
  # Intracluster DE Upregulated genes
  intracluster_active_de_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line, "intracluster_active_de_signatures.rds"))
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  active_signatures <- list()
  
  for(curr_cluster in clusters_of_interest){
    
    cat(curr_cluster, "\n")
    
    inactive_signature <- de_result %>% 
      filter(cluster == curr_cluster & p_val_adj < .05 & avg_log2FC > .5) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    active_signature <- de_result %>% 
      filter(cluster == paste0(curr_cluster,"_resistant") & p_val_adj < .05 & avg_log2FC > .5) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
    
    #Find shared genes between active and inactive componenents
    shared_genes <- intersect(active_signature, inactive_signature)
    
    #Remove shared genes
    active_signature <- active_signature[!active_signature %in% shared_genes]
    
    #Add genes from DE between active and inactive within this cluster
    active_signature <- unique(append(active_signature, intracluster_active_de_signatures[[paste0(curr_cluster,"_active_signature")]]))
    
    active_signatures <- append(active_signatures, list(active_signature))
    
  }
  
  names(active_signatures) <- paste0(curr_cell_line, "_", clusters_of_interest)
  
  all_active_signatures <- append(all_active_signatures, active_signatures)
}

saveRDS(all_active_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/all_active_signatures.rds")
