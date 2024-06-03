args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(ComplexHeatmap)
library(circlize)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

curr_cell_line <- cell_lines[3]

heatmap_list <- list()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  drug_classes <- as.character(unique(data$pathway_level_1))
  drug_classes <- drug_classes[-which(drug_classes == "Vehicle")]
  drug_classes <- drug_classes[-which(drug_classes == "Other")]
  drug_classes <- sort(drug_classes)
  clusters <- as.numeric(levels(data))
  ##############################################################################
  
  pre_cluster_cell_counts <- c()
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    curr_num <- sum(data$pathway_level_1=="Vehicle" & data$Cluster == curr_cluster)
    
    if(curr_num < 10){
      curr_num <- 0
    } 
    
    pre_cluster_cell_counts <- append(pre_cluster_cell_counts, curr_num)
  }
  
  # pre_cluster_cell_counts[pre_cluster_cell_counts < 10] <- 0
  
  total_pre_cells <- sum(data$pathway_level_1=="Vehicle")       
  
  
  heatmap <- matrix(NA, ncol=length(clusters),nrow=length(drug_classes))
  colnames(heatmap) <- paste0("Cluster ", clusters)
  rownames(heatmap) <- drug_classes
  
  curr_drug_class <- drug_classes[3]
  curr_cluster <- 14
  
  for(curr_drug_class in drug_classes){
    cat(curr_drug_class, "\n")
    
    total_curr_cells <- sum(data$pathway_level_1==curr_drug_class)           
    
    for(curr_cluster in clusters){
      
      curr_num <- sum(data$pathway_level_1==curr_drug_class & data$Cluster == curr_cluster)
      
      if(curr_num < 10){
        curr_num <- 0
      } 
      
      contin_table <- matrix(c(curr_num,pre_cluster_cell_counts[curr_cluster],total_curr_cells,total_pre_cells), ncol=2)
      
      if(pre_cluster_cell_counts[curr_cluster] == 0 & curr_num == 0){
        heatmap[curr_drug_class,curr_cluster] <- 1
        
      } else {
        res <- fisher.test(contin_table+1)
        if(res$p.value < 0.05){
          heatmap[curr_drug_class,curr_cluster] <- res$estimate
        } else{
          heatmap[curr_drug_class,curr_cluster] <- res$estimate
        }
      }
    }
  }
  
  
  heatmap <- heatmap[,1:ncol(heatmap) %in% RACs[[curr_cell_line]]]
  
  heatmap <- ifelse(heatmap > 10, 10,heatmap)
  
  col_fun = colorRamp2(c(-1, 0, max(log(heatmap))), c("blue", "white", "red"))
  
  plot_title <- curr_cell_line
  
  ht <- Heatmap(log(heatmap), name="log(OR)", cluster_rows = F,cluster_columns = F,
                column_title = curr_cell_line, col=col_fun, column_names_rot = 45,
                rect_gp = gpar(col = "gray", lwd = .5),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                            labels_gp = gpar(fontsize = 8)))
  
  heatmap_list <- append(heatmap_list, ht)
  
  # 
  
  # 
  # dev.off()
  
}


all_heatmaps <- heatmap_list[[1]] + heatmap_list[[2]] + heatmap_list[[3]]

png(paste0(plotDirectory, "figure_3d.png"),
    width=24, height=10, units= "in", res = 300)

draw(all_heatmaps,column_title="RACs Cell Count Fold Change After Drug Treatment",
     column_title_gp = gpar(fontsize = 30, fontface = "bold"),  padding = unit(c(6, 20, 10, 2), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T,
     ht_gap = unit(1, "cm"))


dev.off()

