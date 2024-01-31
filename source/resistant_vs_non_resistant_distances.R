library(Seurat)
library(tidyverse)

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]


data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")

data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")

data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")


Idents(data) <- data$resistant



DimPlot(data, label = T,split.by = "resistant" ,reduction = "PCA")

resistant_cells <- colnames(data)[data$resistant == "resistant"]

resistant_cells <- sample(resistant_cells, 1000)

resistant_PCs <- data@reductions$PCA@cell.embeddings[resistant_cells,1:10]
resistant_dist <- (dist(resistant_PCs))

non_resistant_cells <- colnames(data)[data$resistant == "nonresistant"]



non_resistant_to_use <- sample(non_resistant_cells, 32825)

non_resistant_PCs <- data@reductions$PCA@cell.embeddings[non_resistant_to_use,1:10]

non_resistant_dist <- (dist(non_resistant_PCs))

wilcox_res <- wilcox.test(resistant_dist, non_resistant_dist)

boxplot(resistant_dist, non_resistant_dist)

if(wilcox_res$p.value < 0.05 & mean(resistant_dist) < mean(non_resistant_dist)){
  count <- count+1
}
# 
# count <- 0
# for(i in 1:100){
#   set.seed(i)
#   cat(i, "\n")
#   
#   non_resistant_to_use <- sample(non_resistant_cells, length(resistant_cells))
#   
#   non_resistant_PCs <- data@reductions$PCA@cell.embeddings[non_resistant_to_use,1:10]
#   
#   non_resistant_dist <- (dist(non_resistant_PCs))
#   
#   wilcox_res <- wilcox.test(resistant_dist, non_resistant_dist)
#   
#   if(wilcox_res$p.value < 0.05 & mean(resistant_dist) < mean(non_resistant_dist)){
#     count <- count+1
#   }
# }
# 
# tests <- count/100
# 
# 
# num_resistant_cells <- length(resistant_cells)
# 
# num_nonresistant_cells <- ncol(data)-num_resistant_cells
# 
# 
# num_nonresistant_cells/num_resistant_cells*1000
# 
# wilcox.test(resistant_dist, non_resistant_dist)
# boxplot(resistant_dist, non_resistant_dist)
# 
# mean(resistant_dist) < mean(non_resistant_dist)
