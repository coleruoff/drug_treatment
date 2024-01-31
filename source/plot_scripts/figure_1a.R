setwd("/data/ruoffcj/projects/drug_treatment/")
library(Seurat)
library(Seurat)
library(ComplexHeatmap)
source("source/cole_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/watermelon_pc9_processed_", genesets_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/watermelon_pc9_processed_", genesets_to_use, "_aucell_thresholds.rds"))
  
  # scores <- scale(scores)
  
  geneset_names <- colnames(scores)
  
  score_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  pvalue_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  
  for(curr_group in groups){
    
    cat(curr_group, "\n")
    
    curr_cell_names <- colnames(data)[data[[ident_to_use]] == curr_group]
    rest_names <- colnames(data)[data[[ident_to_use]] != curr_group]
    
    for(curr_geneset in geneset_names){
      
      curr_scores <- scores[curr_cell_names, curr_geneset]
      
      rest_scores <- scores[rest_names, curr_geneset]
      
      wilcox_res <- wilcox.test(curr_scores,rest_scores)
      
      i <- which(curr_geneset == geneset_names)
      j <- which(curr_group == groups)
      
      score_matrix[i,j] <- mean(curr_scores)
      pvalue_matrix[i,j] <- wilcox_res$p.value
      
    }
  }
  
  scaled_matrix <- t(scale(t(score_matrix)))
  
  colnames(scaled_matrix) <- groups
  rownames(scaled_matrix) <- geneset_names
  
  colnames(score_matrix) <- groups
  rownames(score_matrix) <- geneset_names
  
  return(list(scaled_matrix, pvalue_matrix,score_matrix))
}

watermelon_validation <- function(curr_signature, plot_title=NULL){
  
  if(!exists("watermelon_data")){
    watermelon_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/pc9_processed.rds"))
  }
  
  ################################################################################
  # Calculate FC for each time point
  ################################################################################
  
  Idents(watermelon_data) <- watermelon_data$time_point
  
  time_points <- levels(Idents(watermelon_data))
  
  time_points_fc <- list()
  for(curr_time_point in time_points){
    
    curr_fc_result <- FoldChange(watermelon_data, ident.1 = curr_time_point)
    
    time_points_fc <- append(time_points_fc, list(curr_fc_result))
    
  }
  
  
  names(time_points_fc) <- paste0("pc9_day",time_points,"_fc")
  
  ################################################################################
  # Run GSEA of raj resistance signature for each time point
  ################################################################################
  
  #Create ranked genelists for each time point
  genesets_with_ranks <- list()
  for(curr_fc_res in time_points_fc){
    
    ranks <- curr_fc_res$avg_log2FC
    
    names(ranks) <- rownames(curr_fc_res)
    
    ranks <- sort(ranks, decreasing = T)
    
    genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
    
  }
  
  names(genesets_with_ranks) <- paste0("pc9_day",time_points,"_fc")
  
  
  # Run GSEA
  result <- create_GSEA_matrix(genesets_with_ranks, curr_signature)
  
  gsea_matrix <- result[[1]]
  
  colnames(gsea_matrix) <- paste0("Day ", time_points)
  if(nrow(gsea_matrix) == 1){
    rownames(gsea_matrix) <- ""  
  }
  
  
  if(is.null(plot_title)){
    curr_signature_name <- str_to_title(gsub("_", " ", names(curr_signature)))
    plot_title <- paste0(curr_signature_name, " Enrichment Along Resistance Development in PC9")
  }
  
  
  
  return(gsea_matrix)
}

################################################################################

raj_resistance_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/common_resistance_signature.rds")

raj_watermelon_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_watermelon_resistance_signature.rds")

gsea_matrix <- watermelon_validation(raj_resistance_signature)

plot_title <- "Enrichment of Resistance Signature Along Drug Treatment Time Points"

ht <- Heatmap(gsea_matrix, cluster_columns = F, cluster_rows=F, column_title = plot_title, name="NES",
              column_names_rot = 45, column_title_gp = gpar(fontsize=20, fontface="bold"),column_names_gp = gpar(fontsize=20),
              heatmap_legend_param = list(legend_gp = gpar(fontsize = 20),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                          labels_gp = gpar(fontsize = 12),title_gp = gpar(fontsize = 18, fontface="bold")))


png("/data/ruoffcj/projects/drug_treatment/final_figures/figure_1a.png",width=12,height=5, units = "in", res = 300)

draw(ht)

dev.off()