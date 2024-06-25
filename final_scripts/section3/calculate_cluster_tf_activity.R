# args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
# plotDirectory <- paste0(args[1],"final_figures/")
# setwd(args[1])
library(foreach)
library(doParallel)
library(doMC)
set.seed(42)

registerDoMC(cores=future::availableCores())

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

#################################################################################

cell_lines <- c("A549","K562","MCF7")


RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))
all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

cluster_numbers <- c(19,12,18)
names(cluster_numbers) <- cell_lines

for(curr_cell_line in cell_lines){
  
  data <- all_data[[curr_cell_line]]
  
  # data <- data[,data$pathway_level_1 == "DNA damage & DNA repair"]
  
  num_clusters <- cluster_numbers[[curr_cell_line]]
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_tf_regulons_aucell_scores.rds"))
  #####################################################################################
  
  all_RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))
  
  all_regulons <- colnames(scores)
  
  # tf_heatmap <- matrix(0, ncol=num_clusters,nrow=length(all_regulons))
  
  tf_heatmap <- foreach(curr_regulon = 1:length(all_regulons), .combine = rbind) %dopar% {
    
    cat(curr_regulon, "\n")
    
    curr_row <- c()
    for(curr_cluster in RACs[[curr_cell_line]]){
      
      curr_cluster_names <- colnames(data)[data$Cluster == curr_cluster]
      if(length(curr_cluster_names) == 0){
        curr_value <- 0
      } else{
        
        rest_names <- colnames(data)[!data$Cluster %in% all_RACs[[curr_cell_line]] ]
        
        # Mean score OR for each cluster compared to all non-RACs
        curr_value <- mean(scores[curr_cluster_names,all_regulons[curr_regulon]])/mean(scores[rest_names,all_regulons[curr_regulon]])
      }
      
      
      curr_row <- append(curr_row,curr_value)
      # tf_heatmap[curr_regulon,curr_cluster] <- curr_value
    }
    
    return(curr_row)
  }
  
  
  colnames(tf_heatmap) <- RACs[[curr_cell_line]]
  rownames(tf_heatmap) <- gsub("_regulon", "", all_regulons)
  
  dim(tf_heatmap)
  
  saveRDS(tf_heatmap,paste0(dataDirectory, "processed_data/", curr_cell_line, "_tf_heatmap.rds"))
}


