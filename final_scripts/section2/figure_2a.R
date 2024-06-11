args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))

hallmark_heatmaps <- list()
mp_heatmaps <- list()

curr_cell_line <- cell_lines[1]
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  data <- all_data[[curr_cell_line]]  
  
  hallmark_scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_hallmarks_aucell_scores.rds"))
  mp_scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line, "_processed_filtered_ITH_meta_programs_aucell_scores.rds"))
  
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

  clusters_of_interest <- paste0("Cluster ",RACs[[curr_cell_line]])
  
  rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(curr_hallmark_heatmap) %in% clusters_of_interest,"RAC","Non-RAC")),
                              col = list(RAC = c("RAC" = "orange", "Non-RAC" = "lightblue")), show_annotation_name = F,
                              annotation_legend_param = list(title_gp=gpar(fontsize=22), grid_height=unit(1,"cm"),grid_width=unit(1,"cm"),
                                                             title="Cluster Type", labels_gp = gpar(fontsize = 14)))

    
  
  
  hallmark_ht <- Heatmap(curr_hallmark_heatmap, name="Z-Score", cluster_rows = F, cluster_columns = T,
                bottom_annotation = rac_ha, column_title = curr_cell_line, column_title_side = "top",
                row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=20),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 22),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                            labels_gp = gpar(fontsize = 14)))
  
  
  
  mp_ht <- Heatmap(curr_mp_heatmap, name="Z-Score", cluster_rows = F, cluster_columns = T,
                         bottom_annotation = rac_ha, column_title = curr_cell_line, column_title_side = "top",
                         row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=20),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 22),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                                     labels_gp = gpar(fontsize = 14)))
  
  
  hallmark_heatmaps <- append(hallmark_heatmaps,hallmark_ht)
  mp_heatmaps <- append(mp_heatmaps,mp_ht)
  
}

hallmarks_ht <- hallmark_heatmaps[[1]] + hallmark_heatmaps[[2]] + hallmark_heatmaps[[3]]

png(paste0(plotDirectory,"figure_2a.png"),width=30,height=12, units = "in", res = 300)

draw(hallmarks_ht, column_title = "Cancer Hallmarks Mean AUCell Score", 
     column_title_gp = gpar(fontsize = 36),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T,
     ht_gap = unit(1, "cm"))

dev.off()


png(paste0(plotDirectory,"figure_S2a.png"),
    width=30,height=12, units = "in", res = 300)

mps_ht <- mp_heatmaps[[1]] + mp_heatmaps[[2]] + mp_heatmaps[[3]]

draw(mps_ht, column_title = "Cancer Meta-Programs Mean AUCell Score", 
     column_title_gp = gpar(fontsize = 36),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T,
     ht_gap = unit(1, "cm"))

dev.off()












