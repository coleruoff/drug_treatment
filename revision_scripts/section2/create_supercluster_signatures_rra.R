args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
#################################################################################

cell_lines <- c("A549","K562","MCF7")

num_clusters <- c(19,12,18)
names(num_clusters) <- cell_lines

all_signatures <- list()
for(curr_cell_line in cell_lines){
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  curr_signatures <- list()
  for(curr_cluster in 1:num_clusters[curr_cell_line]){
    
    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de.rds"))
    
    temp <- de_res %>% 
      filter(cluster == curr_cluster & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene)
   
    
    curr_signatures[[curr_cluster]] <- temp
    
  }
  
  all_signatures[[curr_cell_line]] <- curr_signatures
}
################################################################################

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

#Create supercluster signatures
supercluster_signatures <- list()

i <- 3

for(i in 1:length(supercluster_components)){
  
  curr_supercluster_parts <- list()
  #for each cell line
  for(j in names(supercluster_components[[i]])){
    
    #for each component in current cell line
    for(p in supercluster_components[[i]][[j]]){
      
      
      curr_supercluster_parts <- append(curr_supercluster_parts, list(all_signatures[[j]][[p]]))
      
      }
  }
  
  res <- aggregateRanks(curr_supercluster_parts, method = "RRA")
  
  consensus_genes <- subset(res, Score < 0.05)
  
  consensus_genes <- consensus_genes$Name
  
  supercluster_signatures <- append(supercluster_signatures, list(consensus_genes))
}

names(supercluster_signatures) <- paste0("supercluster",1:length(supercluster_signatures),"_signature")

lengths(supercluster_signatures)

saveRDS(supercluster_signatures, paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))
# saveRDS(supercluster_ranks, paste0(dataDirectory, "genesets/supercluster_ranks.rds"))

################################################################################
# Down signatures

all_signatures <- list()
for(curr_cell_line in cell_lines){
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  curr_signatures <- list()
  for(curr_cluster in 1:num_clusters[curr_cell_line]){
    
    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de.rds"))
    
    temp <- de_res %>% 
      filter(cluster == curr_cluster & avg_log2FC < 0 & p_val_adj < 0.05) %>% 
      arrange((avg_log2FC)) %>% 
      pull(gene)
    # select(gene,p_val,avg_log2FC,p_val_adj) 
    
    # ranks <- -log10(temp$p_val_adj+.000000001)*temp$avg_log2FC
    # names(ranks) <- temp$gene
    # 
    # ranks <- sort(ranks,decreasing = T)
    
    
    curr_signatures[[curr_cluster]] <- temp
    
  }
  
  all_signatures[[curr_cell_line]] <- curr_signatures
}
################################################################################

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

#Create supercluster signatures
supercluster_signatures <- list()

i <- 3

for(i in 1:length(supercluster_components)){
  
  curr_supercluster_parts <- list()
  #for each cell line
  for(j in names(supercluster_components[[i]])){
    
    #for each component in current cell line
    for(p in supercluster_components[[i]][[j]]){
      
      
      curr_supercluster_parts <- append(curr_supercluster_parts, list(all_signatures[[j]][[p]]))
      
    }
  }
  
  res <- aggregateRanks(curr_supercluster_parts, method = "RRA")
  
  consensus_genes <- subset(res, Score < 0.05)
  
  consensus_genes <- consensus_genes$Name
  
  supercluster_signatures <- append(supercluster_signatures, list(consensus_genes))
}

names(supercluster_signatures) <- paste0("supercluster",1:length(supercluster_signatures),"_signature")

lengths(supercluster_signatures)

saveRDS(supercluster_signatures, paste0(dataDirectory, "genesets/supercluster_down_signatures.rds"))

