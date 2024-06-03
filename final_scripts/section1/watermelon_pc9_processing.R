args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

library(Seurat)
library(SeuratObject)

cat("SeuratObject: ")
print(packageVersion("SeuratObject"))
cat("\n")

cat("Seurat: ")
print(packageVersion("Seurat"))
cat("\n")


# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

if(!file.exists(paste0(dataDirectory,"processed_data/watermelon_data/"))){
  dir.create(paste0(dataDirectory,"processed_data/watermelon_data/")) 
}

################################################################################

watermelon_pc9.count_matrix <- read.csv(paste0(dataDirectory, "raw_data/watermelon_data/GSE150949_pc9_count_matrix.csv"), row.names = 1)

pc9_metadata <- read.table(paste0(dataDirectory,"raw_data/watermelon_data/GSE150949_metaData_with_lineage.txt"), row.names = 1)

colnames(watermelon_pc9.count_matrix) <- gsub("\\.","-",colnames(watermelon_pc9.count_matrix))

all(rownames(pc9_metadata) == colnames(watermelon_pc9.count_matrix))

watermelon_pc9 <- CreateSeuratObject(counts = watermelon_pc9.count_matrix, meta.data = pc9_metadata, project = "watermelon_pc9_osimertinib")

plot1 <- FeatureScatter(watermelon_pc9, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(watermelon_pc9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

watermelon_pc9 <- NormalizeData(watermelon_pc9)

watermelon_pc9 <- FindVariableFeatures(watermelon_pc9, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(watermelon_pc9)
watermelon_pc9 <- ScaleData(watermelon_pc9, features = all.genes)

# watermelon_pc9 <- RunPCA(watermelon_pc9, features = VariableFeatures(object = watermelon_pc9))

# ElbowPlot(watermelon_pc9)

# watermelon_pc9 <- FindNeighbors(watermelon_pc9, dims = 1:10)

# watermelon_pc9 <- FindClusters(watermelon_pc9, resolution = 0.5)

# watermelon_pc9 <- RunUMAP(watermelon_pc9, dims = 1:10)

# DimPlot(watermelon_pc9, reduction = "umap", group.by = "time_point")

saveRDS(watermelon_pc9, paste0(dataDirectory,"processed_data/watermelon_data/watermelon_pc9_processed.rds"))
# watermelon_pc9 <- readRDS(paste0(dataDirectory,"processed_data/watermelon_data/watermelon_pc9_processed.rds"))

# Run DE for time points
Idents(watermelon_pc9) <- watermelon_pc9$time_point

watermelon_de <- FindAllMarkers(watermelon_pc9)

saveRDS(watermelon_de, paste0(dataDirectory,"de_results/watermelon_time_point_de.rds"))
