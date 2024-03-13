setwd("/data/ruoffcj/projects/drug_treatment/")
source("source/read_in_all_cell_lines.R")
library(ComplexHeatmap)
library(Seurat)
library(circlize)

ident_to_use <- "Cluster"
genesets_to_use <- "hallmarks"

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use){
  
  dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_scores.rds"))
  scores <- scale(scores)
  threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_thresholds.rds"))
  
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
  
  # scaled_matrix <- t(scale(t(score_matrix)))
  scaled_matrix <- scale(score_matrix)
  
  colnames(scaled_matrix) <- groups
  rownames(scaled_matrix) <- geneset_names
  
  colnames(score_matrix) <- groups
  rownames(score_matrix) <- geneset_names
  
  return(list(scaled_matrix, pvalue_matrix,score_matrix))
}

################################################################################

cell_lines <- c("A549","K562","MCF7")

heatmap_matrices <- list()
for(curr_cell_line in cell_lines){

  data <- all_data[[curr_cell_line]]
  
  genesets_name <- "hallmarks"

  ret_mat <- geneset_group_matrix(data, "Cluster", genesets_name)

  #Set all non-significant heatmap cells to 0
  # heatmap_matrix <- ifelse(ret_mat[[2]]<.05,ret_mat[[1]], 0)
  
  heatmap_matrix <- ret_mat[[3]]

  colnames(heatmap_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
  rownames(heatmap_matrix) <- rownames(ret_mat[[1]])

  heatmap_matrices <- append(heatmap_matrices, list(heatmap_matrix))

}

names(heatmap_matrices) <- cell_lines

saveRDS(heatmap_matrices,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_2a_heatmap_matrices.rds")


heatmap_matrices <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_2a_heatmap_matrices.rds")

  
ht_list <- list()
for(curr_cell_line in cell_lines){
  
  heatmap_matrix <- heatmap_matrices[[curr_cell_line]]
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- paste0("Cluster ",RACs[[curr_cell_line]])
  
  rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(heatmap_matrix) %in% clusters_of_interest,"RAC","Non-RAC")),
                              col = list(RAC = c("RAC" = "orange", "Non-RAC" = "lightblue")), show_annotation_name = F,
                              annotation_legend_param = list(title_gp=gpar(fontsize=22), grid_height=unit(1,"cm"),grid_width=unit(1,"cm"),
                                                             title="Cluster Type", labels_gp = gpar(fontsize = 14)))
  
  
  ht <- Heatmap(heatmap_matrix, name="Z-Score", cluster_rows = F, cluster_columns = T,
                bottom_annotation = rac_ha, column_title = curr_cell_line, column_title_side = "top",
                row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=20),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 22),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                            labels_gp = gpar(fontsize = 14)))
  
  
  ht_list <- append(ht_list, ht)
  
}


ht <- ht_list[[1]] + ht_list[[2]] + ht_list[[3]]


# png("/data/ruoffcj/projects/drug_treatment/final_figures/figure_2a.png",width=30,height=12, units = "in", res = 300)

draw(ht, column_title = paste0("Cancer Hallmarks Mean AUCell Score"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T)



dev.off()






