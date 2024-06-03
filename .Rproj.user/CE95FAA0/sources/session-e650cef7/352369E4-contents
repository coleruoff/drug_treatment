library(Seurat)
library(ggplot2)
library(sctransform)
library(tidyverse)

options(future.globals.maxSize = 4000 * 1024^2)

# Load the sample1 dataset
files_to_read <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raw_data/diverse_clonal_fates_data/FM04_BC18_B1/",list.files("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raw_data/diverse_clonal_fates_data/FM04_BC18_B1/"))
sample1.data <- ReadMtx(mtx=files_to_read[3],cells = files_to_read[1], features = files_to_read[2])

files_to_read <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raw_data/diverse_clonal_fates_data/FM04_BC18_B2/",list.files("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raw_data/diverse_clonal_fates_data/FM04_BC18_B2/"))
sample2.data <- ReadMtx(mtx=files_to_read[3],cells = files_to_read[1], features = files_to_read[2])

# Initialize the Seurat object with the raw (non-normalized data).
sample1 <- CreateSeuratObject(counts = sample1.data, project = "10X_Sample1", min.cells = 3, min.features = 200)
sample2 <- CreateSeuratObject(counts = sample2.data, project = "10X_Sample2", min.cells = 3, min.features = 200)

#Calculate percent MT genes
sample1[["percent.mt"]] <- PercentageFeatureSet(object = sample1, pattern = "^MT-")
sample2[["percent.mt"]] <- PercentageFeatureSet(object = sample2, pattern = "^MT-")


#Plot number of features and counts
VlnPlot(object = sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = sample1, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = sample1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))


#Subset cells based on number of features and percent MT
sample1 <- subset(x = sample1, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
sample2 <- subset(x = sample2, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

#Merge two samples together into one Seurat object
s1s2_scTransform <- merge(sample1, y = sample2, add.cell.ids = c("S1", "S2"), project = "S1S2")

Idents(s1s2_scTransform)

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

s1s2_scTransform <- s1s2_scTransform.integrated

s1s2_scTransform <- FindNeighbors(object=s1s2_scTransform, verbose = F)
s1s2_scTransform <- FindClusters(object=s1s2_scTransform, resolution = 0.6, verbose = F)
DimPlot(s1s2_scTransform, reduction = "umap", label = TRUE, pt.size = 0.5)


s1s2_scTransform.markers <- FindAllMarkers(s1s2_scTransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

s1s2_scTransform.markersSubset <- s1s2_scTransform.markers %>% filter(p_val_adj <0.05)


saveRDS(s1s2_scTransform, "//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/FM04_s1s2_integrated.rds")
saveRDS(s1s2_scTransform.markersSubset, "processed_data/diverse_clonal_fates_data/FM04_s1s2_integrated_markers.rds")






genesets <- list()
for(i in unique(s1s2_scTransform.markersSubset$cluster)){
  
  curr_geneset <- s1s2_scTransform.markersSubset %>% 
    filter(cluster == i & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  
  genesets <- append(genesets, list(curr_geneset[1:100]))
}


names(genesets) <- paste0("breast_resistant_cluster", unique(s1s2_scTransform.markersSubset$cluster), "_signature")

saveRDS(genesets, "processed_data/genesets/resistant_breast_signatures_100.rds")





