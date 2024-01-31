library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

#Set cell line
curr_cell_line <- "A549"
#Read in cell line data
data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")

data <- data[,data$Cluster == 9]

                   
resistant <- CreateSeuratObject(counts = data@assays$RNA@counts, project = "resistant_cells", min.cells = 3, min.features = 200)


resistant <- NormalizeData(resistant)

resistant <- FindVariableFeatures(resistant, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(resistant), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(resistant)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(resistant)
resistant <- ScaleData(resistant, features = all.genes)


resistant <- RunPCA(resistant, features = VariableFeatures(object = resistant))

ElbowPlot(resistant)

# devtools::install_version("Matrix",version = "1.6-1.1")

resistant <- FindNeighbors(resistant, dims = 1:10)
resistant <- FindClusters(resistant, resolution = 0.5)

resistant <- RunUMAP(resistant, dims = 1:10)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters


resistant <- AddMetaData(resistant, ifelse(colnames(resistant) %in% active_cell_names, "resistant","nonresistant"), col.name = "resistant")

DimPlot(resistant, reduction = "umap",pt.size = 2)

DimPlot(resistant, group.by = "resistant", label = T, pt.size = 2)


resistant_cell_names <- colnames(resistant)[resistant$resistant == "resistant"]
curr_PCs <- resistant@reductions$pca@cell.embeddings[resistant_cell_names,]

active_dist <- as.matrix(dist(curr_PCs))

# random_cells <- sample(colnames(resistant), length(resistant_cell_names))
random_cells <- colnames(resistant)

curr_PCs <- resistant@reductions$pca@cell.embeddings[random_cells,]

random_dist <- as.matrix(dist(curr_PCs))

active_dist <- as.vector(active_dist)
random_dist <- as.vector(random_dist)

wilcox_res <- wilcox.test(active_dist, random_dist)

boxplot(active_dist, random_dist)


wilcox_res$p.value < 0.05 & median(active_dist) < median(random_dist)
