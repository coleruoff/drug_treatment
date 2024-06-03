library(ComplexHeatmap)
source("source/read_in_all_cell_lines.R")

cell_lines <- c("A549","K562","MCF7")

find_racs <- function(data, curr_cell_line){
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  active_cell_names <- active_cell_names[active_cell_names %in% colnames(data)]
  
  all_clusters <- as.numeric(levels(data))
  
  total_num_active <- length(active_cell_names)
  total_num_inactive <- ncol(data)-total_num_active
  
  ORs <- c()
  pvals <- c()
  
  for(i in all_clusters){
    curr_cluster_names <- colnames(data)[data$Cluster == i]
    
    curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
    curr_num_inactive <- length(curr_cluster_names) - curr_num_active
    
    
    contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
    
    res <- fisher.test(contin_table)
    
    ORs <- append(ORs, res$estimate)
    pvals <- append(pvals, res$p.value)
    
  }
  
  return(all_clusters[ORs > 1.5])
  
}

hallmark_heatmaps <- list()
mp_heatmaps <- list()

curr_cell_line <- cell_lines[1]

for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  data <- all_data[[curr_cell_line]]  
  
  hallmark_scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_hallmarks_aucell_scores.rds"))
  mp_scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_ITH_meta_programs_aucell_scores.rds"))

  # hallmark_scores <- scale(hallmark_scores)
  # mp_scores <- scale(mp_scores)
  
  curr_hallmark_heatmap <- matrix(NA,ncol=0,nrow=ncol(hallmark_scores))
  curr_mp_heatmap <- matrix(NA,ncol=0,nrow=ncol(mp_scores))
  
  drug_classes <- unique(data$pathway_level_1)  

  for(curr_drug_class in drug_classes){
    cat(curr_drug_class, "\n")
    
    # Subset to current drug class
    class_data <- data[,data$pathway_level_1 == curr_drug_class]
    
    # Identify RACs
    curr_racs <- find_racs(class_data, curr_cell_line)
    
    curr_rac_cell_names <- colnames(class_data)[class_data$Cluster %in% curr_racs]
    
    curr_hallmark_scores <- hallmark_scores[rownames(hallmark_scores) %in% curr_rac_cell_names,]
    curr_hallmark_heatmap <- cbind(curr_hallmark_heatmap,colMeans(curr_hallmark_scores))
    
    
    curr_mp_scores <- mp_scores[rownames(mp_scores) %in% curr_rac_cell_names,]
    curr_mp_heatmap <- cbind(curr_mp_heatmap,colMeans(curr_mp_scores))
    
  }
  
  colnames(curr_hallmark_heatmap) <- drug_classes
  colnames(curr_mp_heatmap) <- drug_classes
  
  curr_hallmark_heatmap <- t(scale(t(curr_hallmark_heatmap)))
  curr_mp_heatmap <- t(scale(t(curr_mp_heatmap)))
  
  
  hallmark_ht <- Heatmap(curr_hallmark_heatmap, name="Z-Score", cluster_rows = F, cluster_columns = T,
                column_title = curr_cell_line, column_title_side = "top",
                row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=20),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 22),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                            labels_gp = gpar(fontsize = 14), legend_side="left"))
  
  
  mp_ht <- Heatmap(curr_mp_heatmap, name="Z-Score", cluster_rows = F, cluster_columns = T,
                         column_title = curr_cell_line, column_title_side = "top",
                         row_title = "", column_names_rot = 45,column_title_gp = gpar(fontsize=20),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 22),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                                     labels_gp = gpar(fontsize = 14), legend_side="left"))
    
  
  hallmark_heatmaps <- append(hallmark_heatmaps,hallmark_ht)
  mp_heatmaps <- append(mp_heatmaps,mp_ht)
  
}

hallmarks_ht <- hallmark_heatmaps[[1]] + hallmark_heatmaps[[2]] + hallmark_heatmaps[[3]]

draw(hallmarks_ht, column_title = paste0("Cancer Hallmarks Mean AUCell Score"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T)


mps_ht <- mp_heatmaps[[1]] + mp_heatmaps[[2]] + mp_heatmaps[[3]]

draw(mps_ht, column_title = paste0("Cancer mps Mean AUCell Score"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend=T)













