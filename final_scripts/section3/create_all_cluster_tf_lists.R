args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(ComplexHeatmap)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

#################################################################################s

cell_lines <- c("A549","K562","MCF7")

# all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

cluster_numbers <- c(19,12,18)
names(cluster_numbers) <- cell_lines

all_top_tfs <- list()
all_bottom_tfs <- list()

for(curr_cell_line in cell_lines){
  # data <- all_data[[curr_cell_line]]
  
  num_clusters <- cluster_numbers[curr_cell_line]
  
  tf_heatmap <- readRDS(paste0(dataDirectory, "processed_data/", curr_cell_line, "_tf_heatmap.rds"))
  
  top_tfs <- list()
  bottom_tfs <- list()
  
  for(i in 1:num_clusters){
    
    bottom_thresh <- quantile(tf_heatmap[,i], .5)
    top_thresh <- quantile(tf_heatmap[,i], .5)

    curr_top <- tf_heatmap[tf_heatmap[,i] > top_thresh,i]
    curr_bottom <- tf_heatmap[tf_heatmap[,i] < bottom_thresh,i]
    
    # curr_top <- tf_heatmap[tf_heatmap[,i] > 1.2,i]
    # curr_bottom <- tf_heatmap[tf_heatmap[,i] < .8,i]
    
    # curr_top <- tf_heatmap[tf_heatmap[,i] > 1,i]
    # curr_bottom <- tf_heatmap[tf_heatmap[,i] < 1,i]
    
    # curr_top <- sort(tf_heatmap[,i],decreasing = T)[1:500]
    # curr_bottom <- sort(tf_heatmap[,i],decreasing = F)[1:500]
    
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


