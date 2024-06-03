library(Seurat)
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)


raj_resistant_cancer_type <- "breast"

raj_resistant <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_", raj_resistant_cancer_type, "_processed.rds"))




#Bulkify resistant data
# 
# # Visualize QC metrics as a violin plot
# VlnPlot(raj_resistant, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# 
# plot1 <- FeatureScatter(raj_resistant, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(raj_resistant, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# 
# dim(raj_resistant)
# 
# 
# raj_resistant.filtered <- subset(raj_resistant, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
#          nCount_RNA > 800 &
#          percent.mt < 20)
# 
# dim(raj_resistant.filtered)
DefaultAssay(raj_resistant)
# 
# raj_resistant.filtered <- NormalizeData(raj_resistant.filtered)
aggregated_data <- AggregateExpression(raj_resistant, 
                                       group.by = "seurat_clusters",
                                       assays='RNA',
                                       slot = 'counts',
                                       return.seurat = F)





raj_resistant_bulk <- aggregated_data$RNA

rownames(raj_resistant_bulk) <- rownames(raj_resistant[["RNA"]]@counts)
colnames(raj_resistant_bulk) <- paste0("raj_resistant", 1:ncol(raj_resistant_bulk))



#read in raj control data

if(raj_resistant_cancer_type == "breast"){
  raj_control_cancer_type <- "MDA"
} else if(raj_resistant_cancer_type == "melanoma"){
  raj_control_cancer_type <- "WM989"
}

raj_control <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/memorySeq_data/", raj_control_cancer_type,"_control_counts.rds"))

colnames(raj_control) <- paste0("raj_control", 1:ncol(raj_control))
dim(raj_control)

#Merge resistant and control data

shared_genes <- intersect(rownames(raj_resistant_bulk), rownames(raj_control))

raj_resistant_bulk <- raj_resistant_bulk[shared_genes,]

raj_control <- raj_control[shared_genes,]


all_raj_data <- cbind(raj_resistant_bulk,raj_control)

dim(all_raj_data)

#Create coldata for merged data
coldata <- matrix(c(rep("resistant",ncol(raj_resistant_bulk)), rep("control", ncol(raj_control))), ncol=1)
colnames(coldata) <- c("condition")
rownames(coldata) <- colnames(all_raj_data)
coldata <- data.frame(coldata)

coldata$condition <- factor(coldata$condition)

#Create DESeq Dataset object
dds <- DESeqDataSetFromMatrix(countData = all_raj_data,
                              colData = coldata,
                              design = ~ condition)


#Remove rows with low gene counts
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]


dds$condition <- relevel(dds$condition, ref="control")

dds <- DESeq(dds)

res <- results(dds)
res

resistant_de_genes <- cbind(res@rownames,data.frame(res@listData)) %>% 
  filter(log2FoldChange > 0 & padj < 0.05) %>% 
  arrange(padj) %>% 
  pull("res@rownames")

# res <- as.matrix(res)
saveRDS(res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_", raj_resistant_cancer_type,"_vs_control_deseq_result.rds"))

saveRDS(resistant_de_genes, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type,"_resistant_vs_control_de_genes.rds"))

dds$condition <- relevel(dds$condition, ref="resistant")

dds <- DESeq(dds)

res <- results(dds)
res


control_de_genes <- cbind(res@rownames,data.frame(res@listData)) %>% 
  filter(log2FoldChange > 0 & padj < 0.05) %>% 
  arrange(padj) %>% 
  pull("res@rownames")

# 
# saveRDS(res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_", raj_resistant_cancer_type,"_control_vs_resistant_deseq_result.rds"))
# 
# saveRDS(control_de_genes, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type,"_control_vs_resistant_de_genes.rds"))
# 

###################
#Create genesets for each resistant cluster in raj data

raj_cluster_de_genesets <- list()

for(i in colnames(raj_resistant_bulk)){
  
  cat(i,"\n")
  
  
  cols_to_use <- c(i,colnames(raj_control))
  
  curr_coldata <- coldata %>% 
    filter(rownames(coldata) %in% cols_to_use)
  
  curr_counts <- all_raj_data[,colnames(all_raj_data) %in% cols_to_use]
  
  dds <- DESeqDataSetFromMatrix(countData = curr_counts,
                                colData = curr_coldata,
                                design = ~ condition)
  
  keep <- rowSums(counts(dds)) >= 10
  
  dds <- dds[keep,]
  
  dds$condition <- relevel(dds$condition, ref="control")
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  all_de_genes <- cbind(res@rownames,data.frame(res@listData)) %>% 
    filter(log2FoldChange > 0 & padj < 0.05) %>% 
    arrange(padj) %>% 
    pull("res@rownames")
  
  all_de_genes <- all_de_genes[!grepl("^MT", all_de_genes)]
  all_de_genes <- all_de_genes[!grepl("^RPS", all_de_genes)]
  all_de_genes <- all_de_genes[!grepl("^RPL", all_de_genes)]
  
  raj_cluster_de_genesets <- append(raj_cluster_de_genesets, list(all_de_genes))
}

names(raj_cluster_de_genesets) <- paste0(colnames(raj_resistant_bulk), "_de_genes")


saveRDS(raj_cluster_de_genesets, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_",raj_resistant_cancer_type,"_clusters_vs_control_de_genesets_500.rds"))


###################
#Create genesets for control samples vs each resistant cluster in raj data
# 
# raj_cluster_de_genesets <- list()
# 
# for(i in colnames(raj_resistant_bulk)){
#   
#   cat(i,"\n")
#   
#   
#   cols_to_use <- c(i,colnames(raj_control))
#   
#   curr_coldata <- coldata %>% 
#     filter(rownames(coldata) %in% cols_to_use)
#   
#   curr_counts <- all_raj_data[,colnames(all_raj_data) %in% cols_to_use]
#   
#   dds <- DESeqDataSetFromMatrix(countData = curr_counts,
#                                 colData = curr_coldata,
#                                 design = ~ condition)
#   
#   keep <- rowSums(counts(dds)) >= 10
#   
#   dds <- dds[keep,]
#   
#   dds$condition <- relevel(dds$condition, ref="resistant")
#   
#   dds <- DESeq(dds)
#   
#   res <- results(dds)
#   
#   all_de_genes <- cbind(res@rownames,data.frame(res@listData)) %>% 
#     filter(log2FoldChange > 0 & padj < 0.05) %>% 
#     arrange(padj) %>% 
#     pull("res@rownames")
#   
#   all_de_genes <- all_de_genes[!grepl("^MT", all_de_genes)]
#   all_de_genes <- all_de_genes[!grepl("^RPS", all_de_genes)]
#   all_de_genes <- all_de_genes[!grepl("^RPL", all_de_genes)]
#   
#   raj_cluster_de_genesets <- append(raj_cluster_de_genesets, list(all_de_genes[1:500]))
# }
# 
# names(raj_cluster_de_genesets) <- paste0("control_vs_", colnames(raj_resistant_bulk), "_de_genes")
# 
# 
# saveRDS(raj_cluster_de_genesets, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_control_vs_",raj_resistant_cancer_type,"_clusters_de_genesets_500.rds"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
