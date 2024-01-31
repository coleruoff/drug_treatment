source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")


A549_active <- c(4,9,12,13,14,16,19)

K562_active <- c(5,11)

MCF7_active <-  c(5,7,8,17)

#################################################################################
# Read in Trapnell data
#################################################################################
curr_cell_line <- cell_lines[1]

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

num_clusters <- nlevels(data$Cluster)


cluster_de_results <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/A549_cluster_all_markers_de.rds")

cluster_de_markers <- list()
for(i in unique(de_results$cluster)){
  
  curr_cluster_genes <- cluster_de_results %>% 
    filter(cluster == i & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    pull(gene)
  
  cluster_de_markers <- append(cluster_de_markers, list(curr_cluster_genes))
}


cluster_de_markers[-A549_active]

consensus_resistance <- find_consensus_geneset(cluster_de_markers[-A549_active], 4)








