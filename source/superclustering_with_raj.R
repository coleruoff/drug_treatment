library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#Set cell line
curr_cell_line <- "A549"
#Read in cell line data
A549.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

#Set cell line
curr_cell_line <- "K562"
#Read in cell line data
K562.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

#Set cell line
curr_cell_line <- "MCF7"
#Read in cell line data
MCF7.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))


raj_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds"))


#################################################################################
# Cluster based on mean expression of variable genes
#################################################################################

A549.data <- FindVariableFeatures(A549.data)
A549_variable_genes <- VariableFeatures(A549.data)

K562.data <- FindVariableFeatures(K562.data)
K562_variable_genes <- VariableFeatures(K562.data)

MCF7.data <- FindVariableFeatures(MCF7.data)
MCF7_variable_genes <- VariableFeatures(MCF7.data)

raj_data <- FindVariableFeatures(raj_data, nfeatures = 2000)
raj_variable_genes <- VariableFeatures(raj_data)

all_variable_genes <- unique(c(A549_variable_genes, K562_variable_genes, MCF7_variable_genes,raj_variable_genes))

all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(A549.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(K562.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(MCF7.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(raj_data)]

#################################################################################

cell_lines <- c("A549", "K562", "MCF7", "raj")

heatmap <- matrix(NA, ncol=61, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", 61)


j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  
  if(cell_line != "raj"){
    data <- eval(parse(text=paste0(cell_line, ".data")))
    
    #total_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes,])
    
    clusters <- 1:nlevels(data$Cluster)
    
    for(curr_cluster in clusters){
      cat(curr_cluster,"\n")
      
      cluster_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$Cluster == curr_cluster])
      
      total_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$Cluster != curr_cluster])
      
      #heatmap[,j] <- cluster_mean-total_mean
      heatmap[,j] <- total_mean-cluster_mean
      colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
      
      j <- j + 1
    }
  } else {
    data <- eval(parse(text=paste0(cell_line, "_data")))
    
    #total_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes,])
    
    clusters <- levels(data$seurat_clusters)
    
    for(curr_cluster in clusters){
      cat(curr_cluster,"\n")
      
      cluster_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$seurat_clusters == curr_cluster])
      
      total_mean <- rowMeans(data@assays$RNA@data[all_variable_genes, data$seurat_clusters != curr_cluster])
      
      #heatmap[,j] <- cluster_mean-total_mean
      heatmap[,j] <- total_mean-cluster_mean
      colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
      
      j <- j + 1
    }
  }
}

dim(heatmap)

cor_heatmap <- cor(heatmap, method = "spearman")

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0("A549_",RACs[["A549"]])
clusters_of_interest <- append(clusters_of_interest, paste0("K562_",RACs[["K562"]]))
clusters_of_interest <- append(clusters_of_interest, paste0("MCF7_",RACs[["MCF7"]]))

raj_cols <- colnames(cor_heatmap)[grepl("raj", colnames(cor_heatmap))]

rac_ha <- rowAnnotation(RAC = c(ifelse(colnames(cor_heatmap) %in% clusters_of_interest,"RAC","Non-RAC")),
                        col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue")))


rac_ha <- rowAnnotation(RAC = c(ifelse(colnames(cor_heatmap) %in% clusters_of_interest,"RAC",ifelse(colnames(cor_heatmap) %in% raj_cols, "Raj", "Non-RAC"))),
                        col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue", "Raj"= "pink")))

ht <- Heatmap(cor_heatmap, name="Spearman Correlation", cluster_rows = T, cluster_columns = T,
              right_annotation = rac_ha)


# file_name <- "documents/figures/supercluster_heatmap.png"
# png(file_name, width= 1200, height= 1200, units="px")

draw(ht, column_title="Correlations of Cluster Differential Mean Expression of Most Variable Genes", column_title_gp = gpar(fontsize = 26))

dev.off()






