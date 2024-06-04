args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])
library(Seurat)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(42)

libr_paths <- .libPaths()
cat("LIBRARY DIRECTORY: ", libr_paths[1],"\n")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

for(curr_cancer_type in c("breast","melanoma")){
  if(curr_cancer_type == "breast"){
    sample1.data <- Read10X(data.dir = paste0(dataDirectory, "raw_data/yogesh_et_al/10X_cellRangerOuts/FM04/BC18_B1/filtered_feature_bc_matrix/"))
    sample2.data <- Read10X(data.dir = paste0(dataDirectory, "raw_data/yogesh_et_al/10X_cellRangerOuts/FM04/BC18_B2/filtered_feature_bc_matrix/"))
  } else {
    sample1.data <- Read10X(data.dir = paste0(dataDirectory, "raw_data/yogesh_et_al/10X_cellRangerOuts/FM01/1_1uMPLX/filtered_feature_bc_matrix/"))
    sample2.data <- Read10X(data.dir = paste0(dataDirectory, "raw_data/yogesh_et_al/10X_cellRangerOuts/FM01/2_1uMPLX/filtered_feature_bc_matrix/"))
  }
  
  # Initialize the Seurat object with the raw (non-normalized data).
  sample1 <- CreateSeuratObject(counts = sample1.data, project = "10X_Sample1", min.cells = 3, min.features = 200)
  sample2 <- CreateSeuratObject(counts = sample2.data, project = "10X_Sample2", min.cells = 3, min.features = 200)
  
  sample1[["percent.mt"]] <- PercentageFeatureSet(object = sample1, pattern = "^MT-")
  sample2[["percent.mt"]] <- PercentageFeatureSet(object = sample2, pattern = "^MT-")
  
  #####
  # VlnPlot(object = sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
  # plot1 <- FeatureScatter(object = sample2, feature1 = "nCount_RNA", feature2 = "percent.mt") 
  # plot2 <- FeatureScatter(object = sample2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # CombinePlots(plots = list(plot1,plot2))
  
  if(curr_cancer_type == "breast"){
    sample1 <- subset(x = sample1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 18) 
    sample2 <- subset(x = sample2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 18)
  } else {
    sample1 <- subset(x = sample1, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26) 
    sample2 <- subset(x = sample2, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
  }
  
  s1s2_scTransform <- merge(sample1, y = sample2, add.cell.ids = c("S1", "S2"), project = "S1S2")
  
  s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "orig.ident")
  s1s2_scTransform.list <- s1s2_scTransform.list[c("10X_Sample1", "10X_Sample2")]
  
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
  
  s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
  s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)
  
  # saveRDS(s1s2_scTransform.integrated, file = paste0(dataDirectory,'s1s2ScTransform_50pcs_filter.rds'))
  
  #############################################################################
  s1s2ScTransform <- s1s2_scTransform.integrated
  
  # s1s2ScTransform <- readRDS(paste0(dataDirectory,"s1s2ScTransform_50pcs_filter.rds"))
  s1s2ScTransform <- FindNeighbors(object=s1s2ScTransform, dims=1:50, verbose = FALSE)
  s1s2ScTransform <- FindClusters(object=s1s2ScTransform, resolution = 0.6, verbose = FALSE)
  
  #Combine counts data
  common_genes <- intersect(rownames(s1s2ScTransform@assays$RNA[1]),rownames(s1s2ScTransform@assays$RNA[2]))
  s1s2ScTransform@assays$RNA$counts <- cbind(s1s2ScTransform@assays$RNA[1][common_genes,],s1s2ScTransform@assays$RNA[2][common_genes,])
  
  saveRDS(s1s2ScTransform, paste0(dataDirectory, "processed_data/diverse_clonal_fates_data/raj_resistant_",curr_cancer_type,"_processed.rds"))
}


