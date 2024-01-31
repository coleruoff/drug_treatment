source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

#################################################################################
# Read in Trapnell data
#################################################################################
curr_cell_line <- cell_lines[1]


all_growing_clusters <- list()
for(curr_cell_line in cell_lines){
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  all_clusters <- sort(unique(data$Cluster))
  
  
  total_pre <- sum(data$treatment_stage == 'pre')
  total_post <- sum(data$treatment_stage == 'post')
  
  
  growing_clusters <- c()
  
  fold_changes <- c()
  
  for(curr_cluster in all_clusters){
    
    cluster_pre <- sum(data$Cluster==curr_cluster & data$treatment_stage == 'pre')
    cluster_post <- sum(data$Cluster==curr_cluster & data$treatment_stage == 'post')
    
    
    
    fisher_table <- matrix(c(cluster_post,total_post,cluster_pre,total_pre), nrow = 2,
                           dimnames =
                             list(c("curr cluster", "total"),
                                  c("post", "pre")))
    
    
    fisher_table <- fisher_table+1
    fisher_results <- fisher.test(fisher_table, alternative = "two.sided")
    
    
    fold_changes <- append(fold_changes,log(fisher_results$estimate))
  }
  
  
  growing_clusters <- as.numeric(all_clusters[fold_changes > quantile(fold_changes)[3] & fold_changes > 0])
  all_growing_clusters <- append(all_growing_clusters, list(growing_clusters))
}

names(all_growing_clusters) <- cell_lines

saveRDS(all_growing_clusters, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/growing_clusters.rds")





data@meta.data %>% 
  select(treatment_stage, Cluster) %>% 
  table()


