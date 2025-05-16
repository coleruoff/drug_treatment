source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################

watermelon_pooled.count_matrix <- read.csv(paste0(dataDirectory, "raw_data/oren_data/GSE150949_pooled_watermelon.count.matrix.csv"))

watermelon_pooled.metadata <-  read.csv(paste0(dataDirectory, "raw_data/oren_data/GSE150949_pooled_watermelon.metadata.matrix.csv"), row.names = 1)

colnames(watermelon_pooled.count_matrix) <- gsub("\\.","-",colnames(watermelon_pooled.count_matrix))

all(rownames(watermelon_pooled.metadata) == colnames(watermelon_pooled.count_matrix))

cell_lines <- unique(watermelon_pooled.metadata$cell_line)

saveRDS(cell_lines, paste0(dataDirectory, "oren_cell_lines_names.rds"))

as.data.frame(watermelon_pooled.metadata) %>% 
  dplyr::count(cell_line)

for(curr_cell_line in cell_lines){
  watermelon_pooled <- CreateSeuratObject(counts = watermelon_pooled.count_matrix[,watermelon_pooled.metadata$cell_line == curr_cell_line],
                                          meta.data = watermelon_pooled.metadata[watermelon_pooled.metadata$cell_line == curr_cell_line,],
                                          project = curr_cell_line)
  
  # plot1 <- FeatureScatter(watermelon_pooled, feature1 = "nCount_RNA", feature2 = "percent.mt")
  # plot2 <- FeatureScatter(watermelon_pooled, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # plot1 + plot2
  
  watermelon_pooled <- subset(watermelon_pooled, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  
  watermelon_pooled <- NormalizeData(watermelon_pooled)
  
  watermelon_pooled <- FindVariableFeatures(watermelon_pooled, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(watermelon_pooled)
  watermelon_pooled <- ScaleData(watermelon_pooled, features = all.genes)
  
  watermelon_pooled <- RunPCA(watermelon_pooled, features = VariableFeatures(object = watermelon_pooled))
  
  # ElbowPlot(watermelon_pooled)
  
  watermelon_pooled <- FindNeighbors(watermelon_pooled, dims = 1:10)
  watermelon_pooled <- FindClusters(watermelon_pooled, resolution = 0.5)
  
  watermelon_pooled <- RunUMAP(watermelon_pooled, dims = 1:10)
  
  # DimPlot(watermelon_pooled, reduction = "umap", group.by = "experiment")
  
  file_name <- paste0(dataDirectory,"processed_data/oren_data/", curr_cell_line, "_processed.rds")
  
  saveRDS(watermelon_pooled, file_name)
}
