args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

RACs <- readRDS(paste0(dataDirectory, "processed_data/all_racs.rds"))

hallmark_heatmaps <- list()
mp_heatmaps <- list()


curr_cell_line <- cell_lines[1]
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  data <- all_data[[curr_cell_line]]  
  
  hallmark_scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_hallmarks_aucell_scores.rds"))
  mp_scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_ITH_meta_programs_aucell_scores.rds"))
  
  data <- data[,data$treatment_stage=="post"]
  # hallmark_scores <- t(scale(t(hallmark_scores)))
  # mp_scores <- t(scale(t(mp_scores)))
  
  # Z-score hallmarks
  temp_row_names <- rownames(hallmark_scores)
  temp_col_names <- colnames(hallmark_scores)
  hallmark_scores <- matrix(scale(as.numeric(hallmark_scores)),ncol=50)
  rownames(hallmark_scores) <- temp_row_names
  colnames(hallmark_scores) <- temp_col_names
  
  # Z-score MPs
  temp_row_names <- rownames(mp_scores)
  temp_col_names <- colnames(mp_scores)
  mp_scores <- matrix(scale(as.numeric(mp_scores)),ncol=41)
  rownames(mp_scores) <- temp_row_names
  colnames(mp_scores) <- temp_col_names
  
  curr_hallmark_heatmap <- matrix(NA,ncol=0,nrow=ncol(hallmark_scores))
  curr_mp_heatmap <- matrix(NA,ncol=0,nrow=ncol(mp_scores))
  
  clusters <- sort(unique(data$Cluster))
  
  for(curr_cluster in clusters){
    cat(curr_cluster, "\n")
    
    curr_cluster_cell_names <- colnames(data)[data$Cluster == curr_cluster]
    
    curr_hallmark_scores <- hallmark_scores[rownames(hallmark_scores) %in% curr_cluster_cell_names,]
    curr_hallmark_heatmap <- cbind(curr_hallmark_heatmap,colMeans(curr_hallmark_scores))
    
    
    curr_mp_scores <- mp_scores[rownames(mp_scores) %in% curr_cluster_cell_names,]
    curr_mp_heatmap <- cbind(curr_mp_heatmap,colMeans(curr_mp_scores))
    
  }
  
  colnames(curr_hallmark_heatmap) <- paste0("Cluster ", clusters)
  colnames(curr_mp_heatmap) <- paste0("Cluster ", clusters)
  
  
  curr_hallmark_heatmap <- t(scale(t(curr_hallmark_heatmap)))
  curr_mp_heatmap <- t(scale(t(curr_mp_heatmap)))
  
  
  new_rownames <- sapply(rownames(curr_hallmark_heatmap), FUN = function(x) gsub("HALLMARK_", "", x))
  new_rownames <- sapply(new_rownames, FUN = function(x) gsub("_", " ", x))
  rownames(curr_hallmark_heatmap) <- new_rownames
  
  
  
  hallmark_heatmaps <- append(hallmark_heatmaps, list(curr_hallmark_heatmap))
  mp_heatmaps <- append(mp_heatmaps, list(curr_mp_heatmap))
  
}

hallmark_plots <- list()
mp_plots <- list()
i <- 1
for(i in 1:3){
  
  curr_cell_line <- cell_lines[i]
  
  curr_hallmark_heatmap <- hallmark_heatmaps[[i]]
  curr_mp_heatmap <- mp_heatmaps[[i]]
  
  clusters_of_interest <- paste0("Cluster ",RACs[[curr_cell_line]])
  
  
  # vec1 <- as.vector(cor(curr_hallmark_heatmap[,clusters_of_interest]))
  # 
  # vec2 <- as.vector(cor(curr_hallmark_heatmap[,!colnames(curr_hallmark_heatmap) %in% clusters_of_interest]))
  # 
  # vec1 <- vec1[!vec1 == 1]
  # vec2 <- vec2[!vec2 == 1]
  # 
  # boxplot(vec1,vec2)
  
  
  rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(curr_hallmark_heatmap) %in% clusters_of_interest,"RAC","Non-RAC")),
                              simple_anno_size = unit(1, "mm"),
                              col = list(RAC = c("RAC" = "orange", "Non-RAC" = "lightblue")), show_annotation_name = F,
                              annotation_legend_param = list(title_gp=gpar(fontsize=5), grid_height=unit(2,"mm"),grid_width=unit(2,"mm"),
                                                             title="Cluster Type", labels_gp = gpar(fontsize = 4)))
  
  
  
  hallmark_ht <- Heatmap(curr_hallmark_heatmap, name="Z-Score", cluster_rows = F, cluster_columns = T,
                         bottom_annotation = rac_ha, column_title = curr_cell_line, column_title_side = "top",
                         row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=8),
                         column_names_gp = gpar(fontsize=4),row_names_gp = gpar(fontsize=4),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 5),
                                                     legend_height = unit(2, "mm"), 
                                                     grid_width=unit(2,"mm"),
                                                     labels_gp = gpar(fontsize = 4),
                                                     legend_gp = gpar(lwd = .5)))
  
  
  
  mp_ht <- Heatmap(curr_mp_heatmap, name="Z-Score", cluster_rows = F, cluster_columns = T,
                   bottom_annotation = rac_ha, column_title = curr_cell_line, column_title_side = "top",
                   row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=8),
                   column_names_gp = gpar(fontsize=4),row_names_gp = gpar(fontsize=4),
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 5),
                                               legend_height = unit(2, "mm"), 
                                               grid_width=unit(2,"mm"),
                                               labels_gp = gpar(fontsize = 4),
                                               legend_gp = gpar(lwd = .5)))
  
  
  hallmark_plots <- append(hallmark_plots, list(hallmark_ht))
  mp_plots <- append(mp_plots, list(mp_ht))
  
}

hallmarks_ht <- hallmark_plots[[1]] + hallmark_plots[[2]] + hallmark_plots[[3]]

jpeg(paste0(plotDirectory,"figure_1d.jpg"), width=150, height = 100, units = "mm", res = 600)

draw(hallmarks_ht,
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T,
     ht_gap = unit(2, "mm"))

dev.off()


mps_ht <- mp_plots[[1]] + mp_plots[[2]] + mp_plots[[3]]

jpeg(paste0(plotDirectory,"figure_S1c.jpg"), width=150, height = 100, units = "mm", res = 600)

draw(mps_ht,
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T,
     ht_gap = unit(2, "mm"))

dev.off()








