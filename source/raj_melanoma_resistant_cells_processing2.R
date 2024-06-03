library(Seurat)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(ggrepel)

options(future.globals.maxSize = 4000 * 1024^2)

# Load the sample11umPlx dataset
files_to_read <- paste0(getwd(),"/data/diverse_clonal_fates_data/FM01_A_1uMPLX/",list.files("data/diverse_clonal_fates_data/FM01_A_1uMPLX"))
sample11umPlx.data <- ReadMtx(mtx=files_to_read[3],cells = files_to_read[1], features = files_to_read[2])

files_to_read <- paste0(getwd(),"/data/diverse_clonal_fates_data/FM01_B_1uMPLX/",list.files("data/diverse_clonal_fates_data/FM01_B_1uMPLX"))
sample21umPlx.data <- ReadMtx(mtx=files_to_read[3],cells = files_to_read[1], features = files_to_read[2])

# Initialize the Seurat object with the raw (non-normalized data).
sample11umPlx <- CreateSeuratObject(counts = sample11umPlx.data, project = "10X_Sample1_1uMPLX", min.cells = 3, min.features = 200)
sample21umPlx <- CreateSeuratObject(counts = sample21umPlx.data, project = "10X_Sample2_1uMPLX", min.cells = 3, min.features = 200)

#Calculate percent MT genes
sample11umPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample11umPlx, pattern = "^MT-")
sample21umPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample21umPlx, pattern = "^MT-")


#Plot number of features and counts
VlnPlot(object = sample11umPlx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = sample11umPlx, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = sample11umPlx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))


#Subset cells based on number of features and percent MT
sample11umPlx <- subset(x = sample11umPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
sample21umPlx <- subset(x = sample21umPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

#Merge two samples together into one Seurat object
s1s2_scTransform <- merge(sample11umPlx, y = sample21umPlx, add.cell.ids = c("S1", "S2"), project = "S1S2")


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

#Find integrated data clusters and markers
s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

s1s2_scTransform <- s1s2_scTransform.integrated

s1s2_scTransform <- FindNeighbors(object=s1s2_scTransform, verbose = F)
s1s2_scTransform <- FindClusters(object=s1s2_scTransform, resolution = 0.6, verbose = F)
DimPlot(s1s2_scTransform, reduction = "umap", label = TRUE, pt.size = 0.5)


s1s2_scTransform.markers <- FindAllMarkers(s1s2_scTransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

s1s2_scTransform.markersSubset <- s1s2_scTransform.markers %>% filter(p_val_adj <0.05)


saveRDS(s1s2_scTransform, "processed_data/diverse_clonal_fates_data/FM01_s1s2_integrated.rds")
saveRDS(s1s2_scTransform.markersSubset, "processed_data/diverse_clonal_fates_data/FM01_s1s2_integrated_markers.rds")


# s1_scTransform <- s1s2_scTransform.list[[1]]
# s2_scTransform <- s1s2_scTransform.list[[2]]
# 
#Find sample 1 clusters and markers
# s1_scTransform <- RunPCA(s1_scTransform, verbose = FALSE)
# s1_scTransform <- RunUMAP(s1_scTransform, dims = 1:50)
# 
# s1_scTransform <- FindNeighbors(object=s1_scTransform, verbose = F, dims = 1:50)
# s1_scTransform <- FindClusters(object=s1_scTransform, resolution = 0.6, verbose = F)
# DimPlot(s1_scTransform, reduction = "umap", label = TRUE, pt.size = 0.5)
# 
# 
# s1_scTransform.markers <- FindAllMarkers(s1_scTransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# s1_scTransform.markersSubset <- s1_scTransform.markers %>% filter(p_val_adj <0.05)
# # 
# saveRDS(s1_scTransform, "processed_data/diverse_clonal_fates_data/FM01_s1.rds")
# saveRDS(s1_scTransform.markersSubset, "processed_data/diverse_clonal_fates_data/FM01_s1_markers.rds")

# 
# #Find sample 2 clusters and markers
# s2_scTransform <- RunPCA(s2_scTransform, verbose = FALSE)
# s2_scTransform <- RunUMAP(s2_scTransform, dims = 1:50)
# 
# s2_scTransform <- FindNeighbors(object=s2_scTransform, verbose = F)
# s2_scTransform <- FindClusters(object=s2_scTransform, resolution = 0.6, verbose = F)
# DimPlot(s2_scTransform, reduction = "umap", label = TRUE, pt.size = 0.5)
# 
# 
# s2_scTransform.markers <- FindAllMarkers(s2_scTransform, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# s2_scTransform.markersSubset <- s2_scTransform.markers %>% filter(p_val_adj <0.05)
# 
# saveRDS(s2_scTransform, "processed_data/diverse_clonal_fates_data/FM01_s2.rds")
# saveRDS(s2_scTransform.markersSubset, "processed_data/diverse_clonal_fates_data/FM01_s2_markers.rds")




# write.table(s1s2_scTransform.markers, file=paste0(plotDirectory,'scTransformMarkers_snn06.tsv'), col.names = TRUE, sep='\t')

#s1s2_scTransform.markers = as_tibble(read.table(file = paste0(plotDirectory, "scTransformMarkers_snn06.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

plot = ggplot(data = filter(s1s2_scTransform.markers, !gene %in% list), aes(y=avg_log2FC, x = cluster), color = "gray93", size = 1.5, shape = 16) +
  theme_classic() +
  geom_jitter(width = 0.2, shape = 16)+
  geom_jitter(data = filter(s1s2_scTransform.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC), color = "red", width = 0.2, shape = 16) +
  geom_text_repel(data = filter(s1s2_scTransform.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC, label = gene), color = "red") +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())


plot


list = c("ACTG2", "ACTA2", "MYOCD", "TAGLN", "IFIT2", "OASL", "CCL3", "DDX58", "VCAM1", "PKDCC", "TDO2", "FOXF2", "NGFR", "COL9A3", "S100B", "ITGA6", "GAS7", "MLANA", "SOX10", "MITF", "PMEL", "AXL", "SERPINE1", "BGN")
canonical_resistance_markers <- list(list)
names(list) <- "canonical_resistance_markers"
saveRDS(canonical_resistance_markers, "processed_data/genesets/canonical_resistance_markers.rds")




s1s2_scTransform.markersSubset <- readRDS("processed_data/diverse_clonal_fates_data/FM01_s1s2_integrated_markers.rds")


genesets <- list()
for(i in unique(s1s2_scTransform.markersSubset$cluster)){
  
  curr_geneset <- s1s2_scTransform.markersSubset %>% 
    filter(cluster == i & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  
  genesets <- append(genesets, list(curr_geneset[1:100]))
}


names(genesets) <- paste0("melanoma_resistant_cluster", unique(s1s2_scTransform.markersSubset$cluster), "_signature")

saveRDS(genesets, "processed_data/genesets/melanoma_resistant_signatures_100.rds")

















