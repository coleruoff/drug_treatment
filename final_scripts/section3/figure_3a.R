args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# dataDirectory <- "//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"
#################################################################################

cell_lines <- c("A549", "K562", "MCF7")

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

A549.data <- all_data[["A549"]]
K562.data <- all_data[["K562"]]
MCF7.data <- all_data[["MCF7"]]

A549.data <- A549.data[,A549.data$rac == "rac"]

K562.data <- K562.data[,K562.data$rac == "rac"]

MCF7.data <- MCF7.data[,MCF7.data$rac == "rac"]

#################################################################################
# Cluster based on mean expression of variable genes
#################################################################################

A549.data <- FindVariableFeatures(A549.data)
A549_variable_genes <- VariableFeatures(A549.data)

K562.data <- FindVariableFeatures(K562.data)
K562_variable_genes <- VariableFeatures(K562.data)

MCF7.data <- FindVariableFeatures(MCF7.data)
MCF7_variable_genes <- VariableFeatures(MCF7.data)

all_variable_genes <- unique(c(A549_variable_genes, K562_variable_genes, MCF7_variable_genes))

all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(A549.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(K562.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(MCF7.data)]

length(all_variable_genes)
################################################################################

num_cols <- nlevels(A549.data)+nlevels(K562.data)+nlevels(MCF7.data)

heatmap <- matrix(NA, ncol=num_cols, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", num_cols)

j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  data <- eval(parse(text=paste0(cell_line, ".data")))
  
  clusters <- unique(data$Cluster)
  
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    cluster_mean <- rowMeans(data[["RNA"]]@data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data[["RNA"]]@data[all_variable_genes, data$Cluster != curr_cluster])
    
    heatmap[,j] <- total_mean-cluster_mean
    colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
    
    j <- j + 1
  }
}

cor_heatmap <- cor(heatmap, method = "spearman")


ht <- Heatmap(cor_heatmap, name="Spearman\nCorrelation", cluster_rows = T, cluster_columns = T,
              column_title = "", column_title_side = "bottom",
              row_title_side = "left", row_title_gp = gpar(fontsize=20),
              column_names_rot = 45, 
              row_names_gp = gpar(fontsize=20),
              column_names_gp = gpar(fontsize=18),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 30),legend_height = unit(6, "cm"), grid_width=unit(2,"cm"),
                                          labels_gp = gpar(fontsize = 16)))


png(paste0(plotDirectory, "figure_3a.png"),
    width = 20,height=20, units = 'in',res = 300)

draw(ht, column_title="Correlations of RACs Differential Mean Expression of Most Variable Genes",
     column_title_gp = gpar(fontsize = 30, fontface = "bold"),  padding = unit(c(6, 20, 10, 2), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T)

dev.off()




clustered_heatmap <-  hclust(dist(cor_heatmap, method = "euclidean"), method="complete") 

cluster_groups <- cutree(clustered_heatmap, k=6)

corr_values <- as.vector(cor_heatmap)

corr_values <- corr_values[corr_values<1]

corr_threshold <- quantile(corr_values, probs = 0.90)


supercluster_components <- list()
for(i in unique(cluster_groups)){
  curr_clusters <- names(cluster_groups[cluster_groups == i])
  
  if(sum(cor_heatmap[curr_clusters,curr_clusters] > corr_threshold) == 9){
    
    curr_names <- curr_clusters
    
    curr_clusters <- as.numeric(sapply(curr_clusters, FUN = function(x) gsub("[A-Z|0-9]*_","",x)))
    names(curr_clusters) <- sapply(curr_names, FUN = function(x) gsub("_[0-9]*","",x))
    
    
    supercluster_components <- append(supercluster_components, list(curr_clusters))
  }
}

names(supercluster_components) <- paste0("supercluster", 1:length(supercluster_components))


saveRDS(supercluster_components, paste0(dataDirectory, "processed_data/supercluster_components.rds"))
