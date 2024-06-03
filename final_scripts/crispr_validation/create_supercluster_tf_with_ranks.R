

curr_cell_line <- "A549"
a_heatmap <- readRDS(paste0(dataDirectory, "decoupleR_results/", curr_cell_line, "_tf_heatmap.rds"))

curr_cell_line <- "K562"
k_heatmap <- readRDS(paste0(dataDirectory, "decoupleR_results/", curr_cell_line, "_tf_heatmap.rds"))

curr_cell_line <- "MCF7"
m_heatmap <- readRDS(paste0(dataDirectory, "decoupleR_results/", curr_cell_line, "_tf_heatmap.rds"))


supercluster1_components <- c(9,5,8)
supercluster2_components <- c(14,9,13)

names(supercluster1_components) <- cell_lines
names(supercluster2_components) <- cell_lines

supercluster_components <- list(supercluster1_components,supercluster2_components)

common_genes <- intersect(intersect(rownames(a_heatmap),rownames(k_heatmap)),rownames(m_heatmap))

supercluster_top_tfs <- list()
supercluster_bottom_tfs <- list()

for(i in 1:2){
  sc_heatmap <- cbind(a_heatmap[common_genes,supercluster_components[[i]][1]],
                      k_heatmap[common_genes,supercluster_components[[i]][2]],
                      m_heatmap[common_genes,supercluster_components[[i]][3]])
  
  tf_scores <- sort(rowMeans(sc_heatmap),decreasing = T)
  
  n <- 50
  top_tfs <- names(tf_scores[1:n])
  bottom_tfs <- names(tf_scores[(nrow(sc_heatmap)-(n-1)):nrow(sc_heatmap)])
  
  # bottom_thresh <- quantile(tf_scores, .1)
  # top_thresh <- quantile(tf_scores, .9)
  
  # top_tfs <- names(tf_scores[tf_scores > top_thresh])
  # bottom_tfs <- names(tf_scores[tf_scores < bottom_thresh])
  
  
  supercluster_top_tfs <- append(supercluster_top_tfs,list(top_tfs))
  supercluster_bottom_tfs <- append(supercluster_bottom_tfs,list(bottom_tfs))
}

names(supercluster_top_tfs) <- c("supercluster_1","supercluster_2")
names(supercluster_bottom_tfs) <- c("supercluster_1","supercluster_2")

saveRDS(supercluster_top_tfs, paste0(dataDirectory, "genesets/supercluster_ranks_top_tfs.rds"))
saveRDS(supercluster_bottom_tfs, paste0(dataDirectory, "genesets/supercluster_ranks_bottom_tfs.rds"))



