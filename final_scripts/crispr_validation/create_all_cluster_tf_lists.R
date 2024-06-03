library(ComplexHeatmap)
library(foreach)
library(doParallel)
library(doMC)

registerDoMC(cores=future::availableCores())

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

cell_lines <- c("A549","K562","MCF7")

# all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

non_rac_clusters <- list()
non_rac_clusters[["A549"]] <- list(1,2,3,5,6,7,8,10,11,15,17)
non_rac_clusters[["K562"]] <- list(1,2,3,6,7,8,10,12)
non_rac_clusters[["MCF7"]] <- list(1,2,3,4,6,7,9,10,11,14,15,16,18)

cluster_numbers <- c(19,12,18)
names(cluster_numbers) <- cell_lines

all_top_tfs <- list()
all_bottom_tfs <- list()

for(curr_cell_line in cell_lines){
  # data <- all_data[[curr_cell_line]]
  
  num_clusters <- cluster_numbers[curr_cell_line]
  
  # scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_tf_regulons_aucell_scores.rds"))
  # #####################################################################################
  # 
  # all_regulons <- colnames(scores)
  # 
  # # tf_heatmap <- matrix(0, ncol=num_clusters,nrow=length(all_regulons))
  # 
  # tf_heatmap<- foreach(curr_regulon = 1:length(all_regulons), .combine = rbind) %dopar% {
  # 
  #   cat(curr_regulon, "\n")
  # 
  #   curr_row <- c()
  #   for(curr_cluster in 1:num_clusters){
  # 
  #     curr_cluster_names <- colnames(data)[data$Cluster == curr_cluster]
  #     rest_names <- colnames(data)[data$Cluster %in% non_rac_clusters[[curr_cell_line]]]
  # 
  #     # Mean score OR for each cluster
  #     curr_value <- mean(scores[curr_cluster_names,all_regulons[curr_regulon]])/mean(scores[rest_names,all_regulons[curr_regulon]])
  # 
  #     curr_row <- append(curr_row,curr_value)
  #     # tf_heatmap[curr_regulon,curr_cluster] <- curr_value
  #   }
  # 
  #   return(curr_row)
  # }
  # 
  # 
  # colnames(tf_heatmap) <- 1:num_clusters
  # rownames(tf_heatmap) <- gsub("_regulon", "", all_regulons)
  # 
  # dim(tf_heatmap)
  # 
  # saveRDS(tf_heatmap,paste0(dataDirectory, "decoupleR_results/", curr_cell_line, "_tf_heatmap.rds"))
  tf_heatmap <- readRDS(paste0(dataDirectory, "decoupleR_results/", curr_cell_line, "_tf_heatmap.rds"))
  
  top_tfs <- list()
  bottom_tfs <- list()
  
  for(i in 1:num_clusters){
    
    # bottom_thresh <- quantile(tf_heatmap[,i], .5)
    # top_thresh <- quantile(tf_heatmap[,i], .5)
    # 
    # curr_top <- tf_heatmap[tf_heatmap[,i] > top_thresh,i]
    # curr_bottom <- tf_heatmap[tf_heatmap[,i] < bottom_thresh,i]
    
    # curr_top <- tf_heatmap[tf_heatmap[,i] > 1.2,i]
    # curr_bottom <- tf_heatmap[tf_heatmap[,i] < .8,i]
    
    # curr_top <- tf_heatmap[tf_heatmap[,i] > 1,i]
    # curr_bottom <- tf_heatmap[tf_heatmap[,i] < 1,i]
    
    curr_top <- sort(tf_heatmap[,i],decreasing = T)[1:500]
    curr_bottom <- sort(tf_heatmap[,i],decreasing = F)[1:500]
    
    top_tfs <- append(top_tfs, list(curr_top))
    bottom_tfs <- append(bottom_tfs, list(curr_bottom))
    
    
  }
  
  names(top_tfs) <- 1:num_clusters
  names(bottom_tfs) <- 1:num_clusters
  
  all_top_tfs <- append(all_top_tfs, list(top_tfs))
  all_bottom_tfs <- append(all_bottom_tfs, list(bottom_tfs))
}

names(all_top_tfs) <- cell_lines
names(all_bottom_tfs) <- cell_lines

saveRDS(all_top_tfs, paste0(dataDirectory, "genesets/all_cluster_top_tfs.rds"))
saveRDS(all_bottom_tfs, paste0(dataDirectory, "genesets/all_cluster_bottom_tfs.rds"))


