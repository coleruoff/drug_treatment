setwd("/data/ruoffcj/projects/drug_treatment/")
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

A549.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/A549_processed_filtered2.rds"))
K562.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/K562_processed_filtered2.rds"))
MCF7.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/MCF7_processed_filtered2.rds"))

cell_lines <- c("A549", "K562", "MCF7")
RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

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
##################################################################################
# Type 1 supercluster plots
##################################################################################
#Set cell line
curr_cell_line <- "A549"
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
A549.data <- AddMetaData(A549.data, metadata = ifelse(A549.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
A549.data <- AddMetaData(A549.data, metadata = ifelse(A549.data$rac == "rac" & colnames(A549.data) %in% active_cell_names, "1", ifelse(A549.data$rac == "rac" & (!colnames(A549.data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
A549.data <- AddMetaData(A549.data, metadata = ifelse(A549.data$rac == "rac" & colnames(A549.data) %in% active_cell_names, paste0(A549.data$Cluster, "_1"), ifelse(A549.data$rac == "rac" & (!colnames(A549.data) %in% active_cell_names), paste0(A549.data$Cluster, "_2"), paste0(A549.data$Cluster, "_0"))), col.name = "cell_cluster_group")

A549.data <- A549.data[,A549.data$cell_group == "1"]

#Set cell line
curr_cell_line <- "K562"
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
K562.data <- AddMetaData(K562.data, metadata = ifelse(K562.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
K562.data <- AddMetaData(K562.data, metadata = ifelse(K562.data$rac == "rac" & colnames(K562.data) %in% active_cell_names, "1", ifelse(K562.data$rac == "rac" & (!colnames(K562.data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
K562.data <- AddMetaData(K562.data, metadata = ifelse(K562.data$rac == "rac" & colnames(K562.data) %in% active_cell_names, paste0(K562.data$Cluster, "_1"), ifelse(K562.data$rac == "rac" & (!colnames(K562.data) %in% active_cell_names), paste0(K562.data$Cluster, "_2"), paste0(K562.data$Cluster, "_0"))), col.name = "cell_cluster_group")

K562.data <- K562.data[,K562.data$cell_group == "1"]

#Set cell line
curr_cell_line <- "MCF7"
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(MCF7.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(MCF7.data$rac == "rac" & colnames(MCF7.data) %in% active_cell_names, "1", ifelse(MCF7.data$rac == "rac" & (!colnames(MCF7.data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(MCF7.data$rac == "rac" & colnames(MCF7.data) %in% active_cell_names, paste0(MCF7.data$Cluster, "_1"), ifelse(MCF7.data$rac == "rac" & (!colnames(MCF7.data) %in% active_cell_names), paste0(MCF7.data$Cluster, "_2"), paste0(MCF7.data$Cluster, "_0"))), col.name = "cell_cluster_group")

MCF7.data <- MCF7.data[,MCF7.data$cell_group == "1"]

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

################################################################################

num_cols <- nlevels(A549.data)+nlevels(K562.data)+nlevels(MCF7.data)

heatmap <- matrix(NA, ncol=num_cols, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", num_cols)

j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  data <- eval(parse(text=paste0(cell_line, ".data")))
  
  #total_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes,])
  
  clusters <- unique(data$Cluster)
  
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    cluster_mean <- rowMeans(data@assays$RNA$data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data@assays$RNA$data[all_variable_genes, data$Cluster != curr_cluster])
    
    #heatmap[,j] <- cluster_mean-total_mean
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


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_3c.png"),
    width = 20,height=20, units = 'in',res = 300)

draw(ht, column_title="Correlations of RAC Type 1 Cells Differential Mean Expression of Most Variable Genes",
     column_title_gp = gpar(fontsize = 30, fontface = "bold"),  padding = unit(c(6, 20, 10, 2), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T)

dev.off()