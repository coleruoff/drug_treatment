library(Seurat)

watermelon_pc9.count_matrix <- read.csv("/data/CDSL_hannenhalli/Cole/data/watermelon_data/GSE150949_pc9_count_matrix.csv", row.names = 1)

pc9_metadata <- read.table("/data/CDSL_hannenhalli/Cole/data/watermelon_data/GSE150949_metaData_with_lineage.txt", row.names = 1)

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


watermelon_pc9 <- RunPCA(watermelon_pc9, features = VariableFeatures(object = watermelon_pc9))

ElbowPlot(watermelon_pc9)

watermelon_pc9 <- FindNeighbors(watermelon_pc9, dims = 1:10)
watermelon_pc9 <- FindClusters(watermelon_pc9, resolution = 0.5)

watermelon_pc9 <- RunUMAP(watermelon_pc9, dims = 1:10)

DimPlot(watermelon_pc9, reduction = "umap", group.by = "time_point")

saveRDS(watermelon_pc9, "/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/watermelon_data/watermelon_pc9_processed.rds")
