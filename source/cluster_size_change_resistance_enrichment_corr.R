library(Seurat)
library(tidyverse)
source("source/cole_functions.R")
set.seed(42)

curr_cell_line <- "MCF7"
curr_signature <- raj_early_response_2017


plot_enrichment_vs_size_corr <- function(curr_cell_line, curr_signature){
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  median_cell_count <- data@meta.data %>% 
    count(product_name) %>% 
    pull(n) %>% 
    median()
  
  treatment_stage_counts <- data@meta.data %>% 
    count(treatment_stage) 
  
  total_post <- treatment_stage_counts[1,2]
  total_pre <- treatment_stage_counts[2,2]
  
  
  
  cluster_treatment_stage_counts <- data@meta.data %>% 
    count(Cluster,treatment_stage) 
  
  all_clusters <- sort(as.numeric(unique(data$Cluster)))

  log_ORs <- c()
  for(i in all_clusters){
    curr_cluster_info <- cluster_treatment_stage_counts %>% 
      filter(Cluster == i) 
    
    if(nrow(curr_cluster_info) == 1 | curr_cluster_info[2,3] < (.002*median_cell_count)){
      # log_ORs <- append(log_ORs, log((curr_cluster_info[1,3]/.5)/total_ratio))
      pre_percent <- .5/total_pre
      
      post_percent <- curr_cluster_info[1,3]/total_post
      
      log_ORs <- append(log_ORs, log(post_percent/pre_percent))
      
      
    } else{
      pre_percent <- curr_cluster_info[2,3]/total_pre
      
      post_percent <- curr_cluster_info[1,3]/total_post
      
      log_ORs <- append(log_ORs, log(post_percent/pre_percent))
      
      
    }
  }
  names(log_ORs) <- 1:length(all_clusters)
  
  
  ################################################################################
  trapnell_fc_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_all_clusters_fc.rds"))
  
  genesets_with_ranks <- list()
  for(curr_fc_res in trapnell_fc_results){
    
    ranks <- curr_fc_res$avg_log2FC
    
    names(ranks) <- rownames(curr_fc_res)
    
    ranks <- sort(ranks, decreasing = T)
    
    genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
    
  }
  
  names(genesets_with_ranks) <- paste0(curr_cell_line,"_cluster_", all_clusters)
  
  names(genesets_with_ranks) <- names(trapnell_fc_results)
  
  ################################################################################
  result <- create_GSEA_matrix(genesets_with_ranks, curr_signature)
  
  gsea_results <- result[[2]]
  
  result <- result[[1]]
  
  ################################################################################
  df <- cbind(result[1,],log_ORs)
  colnames(df) <- c("NES","OR")
  rownames(df) <- paste0("Cluster ", 1:nrow(df))
  df <- data.frame(df)
  
  corr_coef <- cor(df$NES,df$OR, method = "spearman")
  corr_coef <- sprintf("%.3f", corr_coef)
  
  curr_signature_name <- str_to_title(gsub("_", " ", names(curr_signature)))
  plt <- ggplot(df,aes(x=NES,y=OR))+
    geom_point()+
    geom_text(label=rownames(df), size=5, vjust=-.5)+
    xlab(paste0(curr_signature_name, " NES"))+
    ylab("Cluster Size Change (log(OR))")+
    ggtitle(paste0(curr_cell_line, " Cluster Size Change vs Cluster Resistance Enrichment (Spearman: ", corr_coef, ")"))+
    theme(axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          title = element_text(size=15))
  
  return(plt)
}
################################################################################


raj_resistant_cancer_type <- "common"

resistant_vs_control_de_genes <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistant_vs_control_de_signature_500.rds"))
resistant_vs_control_de_genes <- list(resistant_vs_control_de_genes)
names(resistant_vs_control_de_genes) <- paste0("raj_", raj_resistant_cancer_type, "_resistant_vs_control_signature") 



p1 <- plot_enrichment_vs_size_corr("A549", raj_early_response_2017)


p1
