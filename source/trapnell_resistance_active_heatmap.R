source("source/cole_functions.R")


create_percentage_heatmap <- function(curr_cell_line, geneset_to_use){
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  num_clusters <- nlevels(data$Cluster)
  
  emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
  
  growing <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/growing_clusters.rds")
  growing <- growing[[curr_cell_line]]
  
  cat(geneset_to_use, "\n")
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
  active_cells <- rownames(scores)[scores>threshold$threshold]
  
  percentage_heatmap <- matrix(NA, nrow=2, ncol=(num_clusters+1))
  
  for(curr_stage in c("pre","post")){
    cat(curr_stage, "\n")
    
    curr_row <- ifelse(curr_stage == "pre", 1, 2)
    
    for(j in 1:num_clusters){
      cat(j, "\n")
      
      curr_cells_names <- rownames(data@meta.data)[data$Cluster == j & data$treatment_stage == curr_stage]
      
      curr_total_cells <- length(curr_cells_names)
      
      if(curr_total_cells < 10){
        
        percentage_heatmap[curr_row,j] <- 0
        
      } else {
        curr_active_cells <- sum(curr_cells_names %in% active_cells)
        
        percentage_heatmap[curr_row,j]  <- curr_active_cells/curr_total_cells
      }
    }
  }
  
  #Add total pre percentage in last column
  pre_cells_names <- rownames(data@meta.data)[data$treatment_stage == "pre"]
  active_pre_cells <- sum(pre_cells_names %in% active_cells)
  
  percentage_heatmap[1, ncol(percentage_heatmap)] <- active_pre_cells/length(pre_cells_names)
  
  #Add total post percentage in last column
  post_cells_names <- rownames(data@meta.data)[data$treatment_stage == "post"]
  active_post_cells <- sum(post_cells_names %in% active_cells)
  
  percentage_heatmap[2, ncol(percentage_heatmap)] <- active_post_cells/length(post_cells_names)
  
  
  colnames(percentage_heatmap)[1:ncol(percentage_heatmap)] <- paste0("Cluster ", 1:ncol(percentage_heatmap))
  colnames(percentage_heatmap)[ncol(percentage_heatmap)] <- "Total"
  rownames(percentage_heatmap) <- c("Pre-Treatment","Post-Treatment")
  
  #################################################################################
  # Plot heatmap
  #################################################################################
  
  plot_title <- paste0("Percentage of ",geneset_to_use, " Active Cells in ",curr_cell_line, " Clusters")
  
  col_fun <-colorRamp2(c(0, 1), c("white", "red"))
  
  growing_ha <- HeatmapAnnotation(Growing = c(ifelse(1:num_clusters %in% growing,"Growing","Non-growing"),"NA"),
                          col = list(Growing = c("Growing" = "red", "Non-growing" = "blue", "NA" = "black")))
  
  ht <- Heatmap(percentage_heatmap, name="Percentage", column_title = plot_title, cluster_rows = F, cluster_columns = F, 
                row_names_side = "left", col = col_fun, column_names_rot = 45,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", percentage_heatmap[i, j]), x, y, gp = gpar(fontsize = 15))},
                bottom_annotation = growing_ha)
  
  return(ht)
}

cell_lines <- c("A549","K562","MCF7")

#################################################################################

geneset_to_use <- "raj_watermelon_resistance_signature"
# geneset_to_use <- "raj_early_response_2017"

#################################################################################

ht <- create_percentage_heatmap("K562","common_resistance_signature")

ht


#################################################################################
#################################################################################
#################################################################################
i <- 1


for(i in 1:length(heatmaps_list)){
  
  temp <- gsub("Percentage of ","", heatmaps_list[[i]]@column_title)
  
  file_name <- paste0(curr_cell_line, "_",gsub(" Active Cells in A549 Clusters", "",temp),"_", raj_resistant_cancer_type,".png",sep = "")
   
  png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/raj_refined_markers_heatmaps/",file_name), width = 2000,height = 1000)
  
  draw(heatmaps_list[[i]])
  
  dev.off()
}
#################################################################################
# Populate tidy data
#################################################################################


df <- data.frame(matrix(NA, nrow=0,ncol=4))

for(curr_geneset in colnames(scores)){
  cat(curr_geneset, "\n")
  
  curr_threshold <- thresholds$threshold[thresholds$gene_set == curr_geneset]
  
  curr_geneset_active_cells <- rownames(scores)[scores[,curr_geneset] > curr_threshold]
  
  
  for(curr_stage in c("pre","post")){
    cat(curr_stage, "\n")
    
    curr_row <- ifelse(curr_stage == "pre", 1, 2)
    
    for(j in 1:num_clusters){
      cat(j, "\n")
      
      curr_cells_names <- rownames(data@meta.data)[data$Cluster == j & data$treatment_stage == curr_stage]
      
      curr_total_cells <- length(curr_cells_names)
      
      if(curr_total_cells < 10){
        
        percentage <- 0
        
      } else {
        curr_active_cells <- sum(curr_cells_names %in% curr_geneset_active_cells)
        
        percentage <- curr_active_cells/curr_total_cells
      }
      
      df <- rbind(df,c(curr_stage,j,curr_geneset,percentage))
      
    }
  }
}

og_df <- df

df <- og_df
colnames(df) <- c("treatment_stage","cluster","geneset","percentage")

df$percentage <- as.numeric(df$percentage)

df$treatment_stage <- factor(df$treatment_stage, levels=c("pre","post"))

df$cluster <- factor(df$cluster,levels = sort(as.numeric(unique(df$cluster))))

ggplot(df, aes(x=cluster,y=percentage, fill=treatment_stage))+
  geom_bar(stat = "identity", position=position_dodge())+
  facet_wrap(~geneset)






