setwd("/data/ruoffcj/projects/drug_treatment/")
library(ComplexHeatmap)
library(Seurat)
library(circlize)
library(tidyverse)
library(ggpubr)

geneset_group_matrix <- function(data, data_name, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/", data_name, "_", genesets_to_use, "_aucell_scores.rds"))
  # threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_thresholds.rds"))
  
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

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
################################################################################
# Set genesets to use as expression space for trapnell data
genesets_name <- "MPs_without_resistance"
raj_genesets_name <- "ITH_meta_programs"
# clean_geneset_name <- "Cancer Hallmarks"
clean_geneset_name <- "ITH Meta-Programs"

cell_lines <- c("A549","K562","MCF7")
curr_cell_line <- "MCF7"

plots <- list()
for(curr_cell_line in cell_lines){
  data <- readRDS(paste0(dataDirectory, "data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  data <- readRDS(paste0(dataDirectory, "data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  ################################################################################
  # Create trapnell geneset heatmap/expression space values
  
  trapnell_mat <- geneset_group_matrix(data, paste0(curr_cell_line,"_processed_filtered"), "cell_cluster_group", genesets_name)
  
  # Heatmap(trapnell_mat[[1]])
  
  ################################################################################
  # Read in raj data
  data2 <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds"))
  
  # Create trapnell geneset heatmap/expression space values
  raj_mat <- geneset_group_matrix(data2, "raj_resistant_breast_processed", "seurat_clusters", raj_genesets_name)
  
  ################################################################################
  
  # trapnell_resistant_mat <- geneset_group_matrix(data, paste0(curr_cell_line,"_processed_filtered"), "resistant_cluster", genesets_name)
  
  
  # Heatmap(trapnell_resistant_mat)
  ################################################################################
  # Z-score heatmaps within each cluster
  
  trapnell_matrix <- scale(trapnell_mat[[3]])
  colnames(trapnell_matrix) <- paste0("trapnell_",colnames(trapnell_matrix))
  # Heatmap(trapnell_matrix)
  
  raj_matrix <- scale(raj_mat[[3]])
  colnames(raj_matrix) <- paste0("raj_",1:ncol(raj_matrix))
  # Heatmap(raj_matrix)
  
  # trapnell_resistant_matrix <- scale(trapnell_resistant_mat[[3]])
  # trapnell_resistant_mat <- trapnell_resistant_matrix[,grepl("resistant",colnames(trapnell_resistant_matrix))]
  # colnames(trapnell_resistant_matrix) <- paste0("trapnell_",colnames(trapnell_resistant_matrix))
  
  ################################################################################
  # Calculate distances between all raj and trapnell clusters
  
  total_matrix <- rbind(t(raj_matrix),t(trapnell_matrix))
  
  Heatmap(total_matrix)
  
  dist_mat <- as.matrix(dist(total_matrix, method="euclidean"))
  subset_dist_mat <- dist_mat[rownames(dist_mat) %in% colnames(raj_matrix),colnames(dist_mat) %in% colnames(trapnell_matrix)]
  
  # Select current RACs 
  # clusters_of_interest <- paste0("trapnell_",RACs[[curr_cell_line]])
  
  # Select RAC resistant subpops
  type_1_clusters <- paste0("trapnell_",RACs[[curr_cell_line]],"_1")
  type_2_clusters <- paste0("trapnell_",RACs[[curr_cell_line]],"_2")
  non_rac_clusters <- paste0("trapnell_",levels(data)[!levels(data) %in% RACs[[curr_cell_line]]],"_0")
  
  # Select distances for RACs and non-RACs
  distances_to_type1 <- as.vector(subset_dist_mat[,type_1_clusters])
  distances_to_type2 <- as.vector(subset_dist_mat[,type_2_clusters])
  non_rac_distances <- as.vector(subset_dist_mat[,non_rac_clusters])
  
  ################################################################################
  # Create dataframe for box plot
  df <- data.frame(cbind(c(distances_to_type1,distances_to_type2,non_rac_distances),c(rep("RAC Type 1", length(distances_to_type1)),rep("RAC Type 2", length(distances_to_type2)),rep("Non-RAC",length(non_rac_distances)))))
  colnames(df) <- c("value","group")
  df$group <- factor(df$group, levels = c("RAC Type 1", "RAC Type 2","Non-RAC"))
  df$value <- as.numeric(df$value)
  
  ################################################################################
  # Plot and save figure as file
  
  # plot_title <- paste0("Cluster Distances to Resistant Clusters in ", clean_geneset_name, " Expression Space (",curr_cell_line,")")
  plot_title <- curr_cell_line
  
  p <- ggboxplot(df, x = "group", y = "value", fill="group",short.panel.labs = FALSE)
  
  my_comparisons <- list(c("RAC Type 1", "RAC Type 2"),c("RAC Type 2", "Non-RAC"),c("RAC Type 1", "Non-RAC"))
  
  p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2,label.y = c(9.7, 10, 10.2), size=8)+
    ggtitle(plot_title)+
    xlab("")+
    ylab("")+
    scale_fill_manual(values=c("red","orange", "lightblue"),name = "Cell Groups")+
    theme(legend.position="right",
          title = element_text(size=20, face = "bold"),
          axis.text = element_text(size=20),
          legend.text = element_text(size=24),
          legend.title = element_text(size=26),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1.5,"cm"))
  
  plots <- append(plots, list(p))
  
  ### Similarity boxplots
  # df$value <- 1/(exp(df$value))
  # 
  # p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  #   ggtitle(plot_title)+
  #   xlab("")+
  #   ylab("Similarity")+
  #   scale_fill_manual(values=c("red","orange", "lightblue"),name = "Cell Groups")+
  #   theme(legend.position="right",
  #         title = element_text(size=20, face = "bold"),
  #         axis.text = element_text(size=30),
  #         legend.text = element_text(size=20))
}


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_1g.png"),
       width=30, height=12, units= "in", res = 300)



figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Distance", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob(paste0("Cell Group Distances to Resistant Clusters in ", clean_geneset_name, " Expression Space"), size=40, face="bold"))


print(p)

dev.off()

