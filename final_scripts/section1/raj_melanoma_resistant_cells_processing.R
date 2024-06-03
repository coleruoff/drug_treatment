library(Seurat)
library(tidyverse)
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- args[1]

if(!file.exists(paste0(dataDirectory,"processed_data/diverse_clonal_fates_data/"))){
  dir.create(paste0(dataDirectory,"processed_data/diverse_clonal_fates_data/")) 
}
################################################################################

# Load the datasets
sample1.data <- Read10X(data.dir = paste0(dataDirectory,"raw_data/yogesh_et_al/10X_cellRangerOuts/FM01/1_1uMPLX/filtered_feature_bc_matrix/"))
sample2.data <- Read10X(data.dir = paste0(dataDirectory,"raw_data/yogesh_et_al/10X_cellRangerOuts/FM01/2_1uMPLX/filtered_feature_bc_matrix/"))

# Initialize the Seurat object with the raw (non-normalized data).
sample1 <- CreateSeuratObject(counts = sample1.data, project = "10X_Sample1_1uMPLX", min.cells = 3, min.features = 200)
sample2 <- CreateSeuratObject(counts = sample2.data, project = "10X_Sample2_1uMPLX", min.cells = 3, min.features = 200)

# Add percent mitochondria metadata
sample1[["percent.mt"]] <- PercentageFeatureSet(object = sample1, pattern = "^MT-")
sample2[["percent.mt"]] <- PercentageFeatureSet(object = sample2, pattern = "^MT-")

# Filter cells and genes
VlnPlot(object = sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = sample1, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = sample2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

sample1 <- subset(x = sample1, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
sample2 <- subset(x = sample2, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

# Integrate two samples
s1s2_scTransform <- merge(sample1, y = sample2, add.cell.ids = c("S1", "S2"), project = "S1S2")

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "orig.ident")
s1s2_scTransform.list <- s1s2_scTransform.list[c("10X_Sample1_1uMPLX", "10X_Sample2_1uMPLX")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]], verbose = FALSE)
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)
s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = FALSE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", 
                                             verbose = FALSE)

# Run dimensionality reduction
s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)

s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

saveRDS(s1s2_scTransform.integrated, file = paste0(dataDirectory,'processed_data/diverse_clonal_fates_data/melanoma_s1s2ScTransform_50pcs_filter.rds'))

################################################################################

s1s2_scTransform <- readRDS(file = paste0(dataDirectory,'processed_data/diverse_clonal_fates_data/melanoma_s1s2ScTransform_50pcs_filter.rds')) 

s1s2_scTransform <- FindNeighbors(object=s1s2_scTransform, verbose = F)

s1s2_scTransform <- FindClusters(object=s1s2_scTransform, resolution = 0.6, verbose = FALSE)

Idents(s1s2_scTransform) <- s1s2_scTransform$seurat_clusters

#Combine counts data
common_genes <- intersect(rownames(s1s2_scTransform@assays$RNA[1]),rownames(s1s2_scTransform@assays$RNA[2]))
s1s2_scTransform@assays$RNA$counts <- cbind(s1s2_scTransform@assays$RNA[1][common_genes,],s1s2_scTransform@assays$RNA[2][common_genes,])

saveRDS(s1s2_scTransform, paste0(dataDirectory, "processed_data/diverse_clonal_fates_data/raj_resistant_melanoma_processed.rds"))

################################################################################






















# rm(sample11umPlx)
# rm(sample11umPlx.data)
# rm(sample21umPlx)
# rm(sample21umPlx.data)
# rm(s1s2_scTransform)
# rm(s1s2_scTransform.list)
# rm(s1s2_scTransform.anchors)
# 
# 
# ###Functionalize Count normalized matrix
# logNormalizedCounts = s1s2_scTransform.integrated[['SCT']]@data #normalized log counts matrix (for counts only, replace "data" by "counts")
# normalizedCounts = s1s2_scTransform.integrated[['SCT']]@counts #normalized counts matrix
# cells_count = colnames(normalizedCounts) #CellIds with Sample number as prefix
# cells_count_cellID = sub("S\\d_", "", cells_count)
# cells_count_sampleNum = gsub("[^S12]", "", cells_count)
# logNormalizedCounts = as_tibble(as.data.frame((t(as.matrix(logNormalizedCounts)))))
# logNormalizedCounts = logNormalizedCounts %>% mutate(cellID = cells_count_cellID,
#                                                      sampleNum = cells_count_sampleNum)
# normalizedCounts = as_tibble(as.data.frame((t(as.matrix(normalizedCounts)))))
# normalizedCounts = normalizedCounts %>% mutate(cellID = cells_count_cellID,
#                                                sampleNum = cells_count_sampleNum)
# 
# write.table(logNormalizedCounts, file=paste0(dataDirectory,'logNormalizedSCTCountsS1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
# write.table(normalizedCounts, file=paste0(dataDirectory,'normalizedSCTCountsS1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
# 
# umapCoordinates = (s1s2_scTransform.integrated[['umap']])@cell.embeddings #UMAP coordinates
# cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
# cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
# cells_UMAP_sampleNum = gsub("[^S12]", "", cells_UMAP)
# umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
# umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
#                                              sampleNum = cells_UMAP_sampleNum)
# write.table(umapCoordinates, file=paste0(dataDirectory,'umapCoordinatesSCT_S1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')