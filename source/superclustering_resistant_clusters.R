library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

#Set cell line
curr_cell_line <- "A549"
#Read in cell line data
A549.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

A549.data <- AddMetaData(A549.data, metadata = ifelse(colnames(A549.data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
A549.data <- AddMetaData(A549.data, metadata = ifelse(A549.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
A549.data <- AddMetaData(A549.data, metadata = ifelse(A549.data$rac == "rac" & A549.data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")

A549.data <- A549.data[,A549.data$resistant == "resistant"]

#Set cell line
curr_cell_line <- "K562"
#Read in cell line data
K562.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

K562.data <- AddMetaData(K562.data, metadata = ifelse(colnames(K562.data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
K562.data <- AddMetaData(K562.data, metadata = ifelse(K562.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
K562.data <- AddMetaData(K562.data, metadata = ifelse(K562.data$rac == "rac" & K562.data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")

K562.data <- K562.data[,K562.data$resistant == "resistant"]
#Set cell line
curr_cell_line <- "MCF7"
#Read in cell line data
MCF7.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(colnames(MCF7.data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(MCF7.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(MCF7.data$rac == "rac" & MCF7.data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")

MCF7.data <- MCF7.data[,MCF7.data$resistant == "resistant"]

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

cell_lines <- c("A549", "K562", "MCF7")

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
    
    cluster_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$Cluster != curr_cluster])
    
    #heatmap[,j] <- cluster_mean-total_mean
    heatmap[,j] <- total_mean-cluster_mean
    colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
    
    j <- j + 1
  }
}

cor_heatmap <- cor(heatmap, method = "spearman")


# clusters_of_interest <- paste0("A549_",RACs[["A549"]], "_resistant")
# clusters_of_interest <- append(clusters_of_interest, paste0("K562_",RACs[["K562"]], "_resistant"))
# clusters_of_interest <- append(clusters_of_interest, paste0("MCF7_",RACs[["MCF7"]], "_resistant"))
# 
# rac_ha <- rowAnnotation(RAC = c(ifelse(colnames(cor_heatmap) %in% clusters_of_interest,"RAC","Non-RAC")),
#                         col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue")))


ht <- Heatmap(cor_heatmap, name="Spearman Correlation", cluster_rows = T, cluster_columns = T)


# file_name <- "documents/figures/supercluster_heatmap.png"
# png(file_name, width= 1200, height= 1200, units="px")

draw(ht, column_title="Correlations of Resistant Cluster Differential Mean Expression of Most Variable Genes", column_title_gp = gpar(fontsize = 26))



dev.off()




