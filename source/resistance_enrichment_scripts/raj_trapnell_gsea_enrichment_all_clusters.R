source("source/cole_functions.R")
library(egg)
library(ggpubr)
require(grid)
set.seed(42)



curr_cell_line <- "A549"

raj_early_response_2017 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_early_response_2017.rds")
raj_resistant_2017 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_resistant_2017.rds")

curr_signature <- raj_early_response_2017

trapnell_gsea_enrichment <- function(curr_cell_line, curr_signature, plot_title=NULL){
  
  emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
  emergent <- emergent[[curr_cell_line]]
  
  non_emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/nonemergent_clusters.rds")
  non_emergent <- non_emergent[[curr_cell_line]]
  
  all_clusters <- sort(c(emergent,non_emergent))
  num_clusters <- length(all_clusters)
  
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
  
  if(is.null(plot_title)){
    curr_signature_name <- str_to_title(gsub("_", " ", names(curr_signature)))
    plot_title <- paste0(curr_signature_name, " Enrichment Along Resistance Development in PC9")
  }
  
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ha <- HeatmapAnnotation(Emergent = ifelse(1:num_clusters %in% emergent,"Emergent","Non-emergent"),
                          col = list(Emergent = c("Emergent" = "red", "Non-emergent" = "blue")))
  
  # png(filename = paste0("output/figures/trapnell_raj_figures/",curr_cell_line,"_all_clusters_raj_GSEA_heatmap.png"),
  # height = 1000, width = 1600, res = 100)
  

  rownames(result) <- ""
  
  ht <- Heatmap((result), cluster_columns = T, column_title = plot_title, name="NES",bottom_annotation = ha, column_names_rot = 45, 
                column_title_gp = gpar(fontsize=20))
  
  
  
  # ht <- Heatmap(t(result), cluster_columns = F, column_title = plot_title, name="NES")
  
  return(list(ht, result))
}

################################################################################

temp <- list("test"=early_drug_markers[1:200])

ht <- trapnell_gsea_enrichment("A549", temp)


ht[[1]]

wilcox.test(ht[[2]][14:19],ht[[2]][1:14])
boxplot(ht[[2]][14:19],ht[[2]][1:14])
