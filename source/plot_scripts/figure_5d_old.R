setwd("/data/ruoffcj/projects/drug_treatment/")
library(ComplexHeatmap)
library(Seurat)
library(circlize)
library(ggpubr)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"

# ident_to_use <- "cell_cluster_group"
# genesets_to_use <- "yeast_upregulated_orthologs"

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_thresholds.rds"))
  
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

################################################################################

curr_cell_line <- "MCF7"

cell_lines <- c("A549","K562","MCF7")

plots <- list()

heatmap_matrices <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line, "\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  
  genesets_name <- "ecoli_AMR_genesets_orthologs"
  
  ret_mat <- geneset_group_matrix(data, "cell_cluster_group", genesets_name)
  
  #Set all non-significant heatmap cells to 0
  heatmap_matrix <- ifelse(ret_mat[[2]]<.05,ret_mat[[1]], 0)
  
  colnames(heatmap_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
  rownames(heatmap_matrix) <- rownames(ret_mat[[1]])
  
  heatmap_matrix_up <- heatmap_matrix[!grepl("down",rownames(heatmap_matrix)),]
  
  rownames(heatmap_matrix_up) <- sapply(X = rownames(heatmap_matrix_up), FUN = function(x) {strsplit(x,split = "_")[[1]][1]})
  
  heatmap_matrices <- append(heatmap_matrices, list(heatmap_matrix_up))
  
  df <- data.frame(cbind(colnames(heatmap_matrix_up),heatmap_matrix_up[1,]))
  
  colnames(df) <- c("cluster","value")
  df$value <- as.numeric(df$value)
  
  rac_ha <- HeatmapAnnotation(cell_group = c(ifelse(grepl("_0",colnames(heatmap_matrix_up)),"Non-RAC", ifelse(grepl("_1",colnames(heatmap_matrix_up)), "RAC Type 1","RAC Type 2"))),
                              col = list(cell_group = c("RAC Type 1" = "red","RAC Type 2" = "orange", "Non-RAC" = "lightblue")), show_annotation_name = F,
                              annotation_legend_param = list(legend_gp=gpar(fontsize=20), grid_height=unit(1,"cm"),grid_width=unit(1,"cm"),
                                                             title="Cluster Type", labels_gp = gpar(fontsize = 14)))
  
  colnames(heatmap_matrix_up) <- gsub("_0", "", colnames(heatmap_matrix_up))
  colnames(heatmap_matrix_up) <- gsub("_1", " (Type 1)", colnames(heatmap_matrix_up))
  colnames(heatmap_matrix_up) <- gsub("_2", " (Type 2)", colnames(heatmap_matrix_up))
  
  
  
  
  
  ht <- Heatmap(heatmap_matrix_up, name="Z-Score", cluster_rows = F, cluster_columns = T,
                bottom_annotation = rac_ha, column_title = "", column_title_side = "bottom",
                row_title = "Antimicrobial Drugs\n", row_title_side = "left", row_title_gp = gpar(fontsize=25),
                row_names_side = "left", column_names_rot = 45, 
                row_names_gp = gpar(fontsize=20),
                column_names_gp = gpar(fontsize=20),
                heatmap_legend_param = list(legend_gp = gpar(fontsize = 12),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                            labels_gp = gpar(fontsize = 14)))
  
  
  p <- grid.grabExpr(draw(ht, column_title = curr_cell_line, 
       column_title_gp = gpar(fontsize = 25, fontface = "bold"),  padding = unit(c(2, 2, 10, 2), "mm"),
       heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T))
  

  plots <- append(plots,list(p))  
}


# saveRDS(heatmap_matrices, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_5d_heatmap_matrices.rds")
heatmap_matrices <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/figure_data/figure_5d_heatmap_matrices.rds")

cell_lines <- c("A549","K562","MCF7")

plots <- list()
for(i in 1:length(heatmap_matrices)){
  
  heatmap_matrix_up <- heatmap_matrices[[1]]
  
  curr_cell_line <- cell_lines[i]
  
  df <- data.frame(cbind(colnames(heatmap_matrix_up),heatmap_matrix_up[1,]))
  
  colnames(df) <- c("cluster","value")
  df$value <- as.numeric(df$value)
  
  rac_ha <- HeatmapAnnotation(cell_group = c(ifelse(grepl("_0",colnames(heatmap_matrix_up)),"Non-RAC", ifelse(grepl("_1",colnames(heatmap_matrix_up)), "RAC Type 1","RAC Type 2"))),
                              col = list(cell_group = c("RAC Type 1" = "red","RAC Type 2" = "orange", "Non-RAC" = "lightblue")), show_annotation_name = F, show_legend = F)
  
  # colnames(heatmap_matrix_up) <- gsub("_0", "", colnames(heatmap_matrix_up))
  # colnames(heatmap_matrix_up) <- gsub("_1", " (Type 1)", colnames(heatmap_matrix_up))
  # colnames(heatmap_matrix_up) <- gsub("_2", " (Type 2)", colnames(heatmap_matrix_up))
  
  colnames(heatmap_matrix_up) <- gsub("_0", "", colnames(heatmap_matrix_up))
  colnames(heatmap_matrix_up) <- gsub("_1", "", colnames(heatmap_matrix_up))
  colnames(heatmap_matrix_up) <- gsub("_2", "", colnames(heatmap_matrix_up))
  
  
  
  ht <- Heatmap(heatmap_matrix_up, name="Z-Score", cluster_rows = F, cluster_columns = T,
                bottom_annotation = rac_ha, column_title = "", column_title_side = "bottom",
                row_title = "", row_title_side = "left", row_title_gp = gpar(fontsize=25),
                row_names_side = "left", column_names_rot = 45, 
                row_names_gp = gpar(fontsize=20),
                column_names_gp = gpar(fontsize=20), show_heatmap_legend = F)
 
 

  
  p <- grid.grabExpr(draw(ht, column_title = curr_cell_line, 
                          column_title_gp = gpar(fontsize = 25),  padding = unit(c(2, 2, 10, 2), "mm"),
                          heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T))
  
  
  plots <- append(plots,list(p))  
}

# Create legends
anno_legend <- Legend(labels = c("RAC Type 1","RAC Type 2","Non-RAC"), title = "Cluster Type",
                      legend_gp = gpar(fill=c("red","orange","lightblue",fontsize=25)),
                      grid_height=unit(1,"cm"),grid_width=unit(1,"cm"), 
                      labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize=30))

col_fun = colorRamp2(c(-2, 0.5, 2), c("blue", "white", "red"))
heatmap_legend <- Legend(col_fun=col_fun,title="Z-Score",legend_gp = gpar(fontsize = 12),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                         labels_gp = gpar(fontsize = 14), title_gp =gpar(fontsize=30))

legends <- packLegend(heatmap_legend,anno_legend)

legends <- grid.grabExpr(draw(legends))

# Arrange plots and legends
plots <- append(plots, list(legends))
lay <- rbind(c(1,1,1,1,4),
             c(2,2,2,2,4),
             c(3,3,3,3,4))





png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_5d.png"),
    width=20, height=25, units="in",res = 300)

grid.arrange(grobs = plots, layout_matrix = lay,top=textGrob("E. Coli Antimicrobial Resistance Orthologs Mean Cluster AUCell Score",gp=gpar(fontsize=40, fontface="bold")),
             left=textGrob("Antimicrobial Drugs\n", rot=90,gp=gpar(fontsize=30, fontface="bold")))

dev.off()





