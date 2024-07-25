args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(ComplexHeatmap)
library(circlize)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

curr_cell_line <- cell_lines[3]

heatmap_list <- list()

curr_cell_line <- "A549"

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  drug_classes <- as.character(unique(data$pathway_level_1))
  drug_classes <- drug_classes[-which(drug_classes == "Vehicle")]
  drug_classes <- drug_classes[-which(drug_classes == "Other")]
  drug_classes <- sort(drug_classes)
  
  clusters <- c(supercluster_components[[1]][[curr_cell_line]],supercluster_components[[2]][[curr_cell_line]])
  ##############################################################################
   
  heatmap <- matrix(NA, ncol=length(clusters),nrow=length(drug_classes))
  colnames(heatmap) <- paste0("Cluster ", clusters)
  rownames(heatmap) <- drug_classes
  
  total_cells <- ncol(data) 
  
  curr_cell_line <- "A549"
  curr_cluster <- 14
  curr_drug_class <- drug_classes[1]
  
  for(curr_drug_class in drug_classes){
    cat(curr_drug_class, "\n")
    
    total_drug_cells <- sum(data$pathway_level_1==curr_drug_class)  
     
    for(curr_cluster in clusters){
      
      curr_num <- sum(data$pathway_level_1==curr_drug_class & data$Cluster == curr_cluster)
      cluster_total <- sum(data$Cluster == curr_cluster)
      
      i <- which(curr_cluster == clusters)
      heatmap[curr_drug_class,i] <- (curr_num/cluster_total)/(total_drug_cells/total_cells)
      
    }
  }
  
  heatmap_list <- append(heatmap_list, list(heatmap))
  
}


sc_heatmap_list <- list()
for(i in 1:2){
  curr_sc_heatmap <- cbind(heatmap_list[[1]][,i],heatmap_list[[2]][,i],heatmap_list[[3]][,i])
  
  curr_sc_heatmap[curr_sc_heatmap == 0] <- NA
  
  colnames(curr_sc_heatmap) <- cell_lines
  
  heatmap_title <- paste0("Supercluster ", i)
  
  ht <- Heatmap(log(curr_sc_heatmap), name="log(OR)", cluster_rows = F,cluster_columns = F,
                              column_title = heatmap_title,  column_names_rot = 45,column_names_gp = gpar(fontsize= 6),
                              rect_gp = gpar(col = "gray", lwd = .5),column_title_gp = gpar(fontsize=6),
                              row_names_gp = gpar(fontsize = 4),
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 8),legend_height = unit(5, "mm"), grid_width=unit(1,"mm"),
                                                          labels_gp = gpar(fontsize = 6)))
  
  sc_heatmap_list <- append(sc_heatmap_list, list(ht))
}



all_heatmaps <- sc_heatmap_list[[1]] + sc_heatmap_list[[2]]

tiff(paste0(plotDirectory,"figure_2e.tiff"), width=90, height = 40, units = "mm", res = 1000)

draw(all_heatmaps,
     heatmap_legend_side = "left",merge_legend=T,
     ht_gap = unit(1, "mm"))

dev.off()

# png(paste0(plotDirectory, "figure_2e.png"),
#     width=20, height=10, units= "in", res = 300)
# 
# draw(all_heatmaps,padding = unit(c(6, 20, 10, 100), "mm"),
#      heatmap_legend_side = "left",merge_legend=T,
#      ht_gap = unit(1, "cm"))
# 
# 
# dev.off()

