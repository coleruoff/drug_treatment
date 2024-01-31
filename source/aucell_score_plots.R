library(ggpubr)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


# curr_cell_line <- "A549"
# ident_to_use <- "Cluster"
# geneset_to_use <- "raj_watermelon_resistance_signature"


create_aucell_plots <- function(data, ident_to_use, geneset_to_use){
  
  all_boxplots <- list()
  all_heatmaps <- list()
  
  # num_clusters <- nlevels(data$Cluster)
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  # emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
  # emergent_clusters <- emergent[[curr_cell_line]]
  # data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% emergent_clusters, data$Cluster,"nonemergent"), col.name = "emergent")
  # growing <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/growing_clusters.rds")
  # growing_clusters <- growing[[curr_cell_line]]
  # data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% growing_clusters, data$Cluster,"nongrowing"), col.name = "growing")
  
  cat(geneset_to_use, "\n")
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
  
  # Create list of lists containing active cells for each geneset
  all_active_cells <- list()
  
  curr_geneset <- colnames(scores)[1]
  for(curr_geneset in colnames(scores)){
    
    cat(curr_geneset, "\n")
    active_cells <- rownames(scores)[scores[,curr_geneset] > threshold$threshold[threshold$gene_set == curr_geneset]]
    
    # all_active_cells <- append(all_active_cells, list(curr_active_cells))
    
    data <- AddMetaData(data, metadata = scores[,curr_geneset], col.name = "curr_geneset")
    
    # df <- data@meta.data %>% 
    #   dplyr::select(ident_to_use, curr_geneset, emergent, dose, treatment_stage) %>% 
    #   group_by(emergent) %>% 
    #   mutate("emergent_median_score" = median(curr_geneset)) %>% 
    #   group_by(ident_to_use) %>% 
    #   mutate("group_median_score" = median(curr_geneset))
    
    df <- data@meta.data %>% 
      dplyr::select(ident_to_use, curr_geneset, dose, treatment_stage) %>% 
      group_by(eval(parse(text=ident_to_use))) %>% 
      mutate("group_median_score" = median(curr_geneset))
    
    
    
    p <- ggboxplot(df, x = ident_to_use, y = "curr_geneset", fill="group_median_score")+
      scale_fill_gradient(low="white", high="red")+
      ggtitle(paste0(curr_cell_line, " Cluster Scores for ", curr_geneset))
    #  Add p-value
    p2 <- p + stat_compare_means()
    
    all_boxplots <- append(all_boxplots, list(p2))
    
    
    percentage_heatmap <- matrix(NA, nrow=2, ncol=(num_groups+1))
    
    for(curr_stage in c("pre","post")){
      cat(curr_stage, "\n")
      
      curr_row <- ifelse(curr_stage == "pre", 1, 2)
      
      for(j in 1:num_groups){
        cat(j, "\n")
        
        curr_cells_names <- rownames(data@meta.data)[data[[ident_to_use]] == j & data$treatment_stage == curr_stage]
        
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
    
    plot_title <- paste0("Percentage of ",curr_geneset, " Active Cells in ",curr_cell_line, " Clusters")
    
    col_fun <-colorRamp2(c(0, 1), c("white", "red"))
    
    # growing_ha <- HeatmapAnnotation(Growing = c(ifelse(1:num_groups %in% growing_clusters,"Growing","Non-growing"),"NA"),
    #                                 col = list(Growing = c("Growing" = "red", "Non-growing" = "blue", "NA" = "black")))
    
    ht <- Heatmap(percentage_heatmap, name="Percentage", column_title = plot_title, cluster_rows = F, cluster_columns = F, 
                  row_names_side = "left", col = col_fun, column_names_rot = 45,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", percentage_heatmap[i, j]), x, y, gp = gpar(fontsize = 15))})
    
    all_heatmaps<- append(all_heatmaps, ht)
  }
  
  names(all_heatmaps) <- colnames(scores)
  names(all_boxplots) <- colnames(scores)
  
  
  return(list("heatmaps" = all_heatmaps, "boxplots" = all_boxplots))
}

###################################################################################

curr_cell_line <- "K562"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")

results <- create_aucell_plots(data, "Cluster", "raj_watermelon_resistance_signature")


results$boxplots$raj_watermelon_resistance_signature
results$heatmaps$raj_watermelon_resistance_signature


#################################################################################

# finding RACs
percentage_mat <- results$heatmaps$raj_watermelon_resistance_signature@matrix



threshold <- .089
plot(density(percentage_mat[2,]),
     xlab="Percent Active Cells", #Change the x-axis label
     ylab="Density", #y-axis label
     main=curr_cell_line)
abline(v=threshold, col= "red",lwd = 4)


colnames(percentage_mat)[percentage_mat[2,]>=threshold]


###################################################################################




























curr_cell_line <- cell_lines[1]

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

geneset_to_use <- "raj_watermelon_resistance_signature"
# geneset_to_use <- "raj_early_response_2017"
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
active_cells <- rownames(scores)[scores>threshold$threshold]

emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
emergent_clusters <- emergent[[curr_cell_line]]
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% emergent_clusters, data$Cluster,"nonemergent"), col.name = "emergent")

data <- AddMetaData(data, metadata = scores[,1], col.name = "raj_resistant")

data <- data[,data$treatment_stage == "post"]


df <- data@meta.data %>% 
  select(Cluster, raj_resistant, emergent, dose, treatment_stage) %>% 
  group_by(emergent) %>% 
  mutate("emergent_median_score" = median(raj_resistant)) %>% 
  group_by(Cluster) %>% 
  mutate("cluster_median_score" = median(raj_resistant))



p <- ggboxplot(df, x = "Cluster", y = "raj_resistant", fill="cluster_median_score")+
  scale_fill_gradient(low="white", high="red")+
  ggtitle(paste0(curr_cell_line, " Cluster Scores for ", geneset_to_use))
#  Add p-value
p + stat_compare_means()


cluster_active_percentages <- c()
for(i in sort(unique(df$Cluster))){
 
  cluster_active_percentages <- append(cluster_active_percentages,sum(data$Cluster == i & colnames(data) %in% active_cells)/sum(data$Cluster == i) )
  
}


boxplot(cluster_active_percentages[1:13],cluster_active_percentages[14:19])

boxplot(cluster_active_percentages[-9],cluster_active_percentages[9])

plot(density(cluster_active_percentages))
abline(v=.05)
(1:length(unique(df$Cluster)))[cluster_active_percentages > .05]



percent_df <- data.frame(cbind(1:length(cluster_active_percentages),cluster_active_percentages))

percent_df$cluster <- factor(as.character(percent_df$cluster), levels = as.character(percent_df$cluster))

colnames(percent_df) <- c("cluster","value")

ggplot(percent_df)+
  geom_col(aes(x=cluster,y=value))+
  xlab("Cluster")+
  ylab("Percentage of Active Cells")+
  ggtitle(paste0("Percentage of Active Cells in Post-Treatment ", curr_cell_line))+
  theme(axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        title = element_text(size=25))


for(i in sort(unique(df$Cluster))){
  pre_scores <- df %>% 
    filter(Cluster == i & treatment_stage=="pre") %>% 
    pull(raj_resistant)
  
  post_scores <- df %>% 
    filter(Cluster == i & treatment_stage=="post") %>% 
    pull(raj_resistant)
  
  boxplot(pre_scores,post_scores)
  
  if(length(pre_scores) > 10){
    
    res <- wilcox.test(pre_scores,post_scores, alternative = "less")
    
    if(res$p.value < 0.05){
      cat(i,"\n")
    }
  }
}





p <- ggboxplot(df, x = "emergent", y = "raj_resistant", fill="emergent_median_score")+
  scale_fill_gradient(low="white", high="red")+
  ggtitle(paste0(curr_cell_line, " Cluster Scores for ", geneset_to_use))
#  Add p-value
p + stat_compare_means()




p <- ggboxplot(df, x = "active_cluster", y = "raj_resistant", fill="active_cluster_median_score")+
  scale_fill_gradientn(colours=c("blue","white","red"),na.value = "transparent",breaks=c(-1.5,0,1.5))+
  ggtitle(paste0(curr_cell_line, " Cluster Scores for ", geneset_to_use))
#  Add p-value
p + stat_compare_means()




cluster_median_scores <- data@meta.data %>% 
  select(Cluster, raj_resistant, emergent) %>% 
  group_by(Cluster) %>% 
  mutate("cluster_median_score" = median(raj_resistant)) %>% 
  arrange(Cluster) %>% 
  select(cluster_median_score) %>% 
  unique()

