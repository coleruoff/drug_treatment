library(Seurat)
library(tidyverse)
library(DESeq2)
source("source/cole_functions.R")
set.seed(42)

################################################################################
#Read in raj resistant data
################################################################################
raj_resistant_cancer_type <- "melanoma"

raj_resistant <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_", raj_resistant_cancer_type, "_processed.rds"))

################################################################################
#Read in raj control data
################################################################################
if(raj_resistant_cancer_type == "breast"){
  raj_control_cancer_type <- "MDA"
} else if(raj_resistant_cancer_type == "melanoma"){
  raj_control_cancer_type <- "WM989"
}

raj_control <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/memorySeq_data/", raj_control_cancer_type,"_control_counts.rds"))

colnames(raj_control) <- paste0("raj_control", 1:ncol(raj_control))

################################################################################
# Find resistant cluster markers
################################################################################

# cluster_de_results <- FindAllMarkers(raj_resistant)

# saveRDS(cluster_de_results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_", raj_resistant_cancer_type, "_cluster_de.rds"))

# cluster_de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_", raj_resistant_cancer_type, "_cluster_de.rds"))

# cluster_de_markers <- list()
# for(i in unique(cluster_de_results$cluster)){
#   curr_cluster_genes <- cluster_de_results %>% 
#     filter(cluster == i & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     pull(gene)
#   
#   cluster_de_markers <- append(cluster_de_markers, list(curr_cluster_genes))
# }



# 
# cluster_de_consensus <- find_consensus_geneset(cluster_de_markers, 3)
# 
# cluster_de_consensus <- sapply(cluster_de_consensus, FUN = function(x){  gsub("\\.*", "",x)   })
# 
# names(cluster_de_consensus) <- NULL

################################################################################
# Create resistant vs control DE results
################################################################################
aggregated_data <- AggregateExpression(raj_resistant, 
                                       group.by = "seurat_clusters",
                                       assays='RNA',
                                       slot = 'counts',
                                       return.seurat = F)

raj_resistant_bulk <- aggregated_data$RNA

rownames(raj_resistant_bulk) <- rownames(raj_resistant[["RNA"]]@counts)
colnames(raj_resistant_bulk) <- paste0("raj_resistant", 1:ncol(raj_resistant_bulk))


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

saveRDS(res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_", raj_resistant_cancer_type, "_resistant_vs_control_de_res.rds"))

################################################################################
# Create cancer type specific resistance signatures
################################################################################
signature_length <- 200

resistant_vs_control_de_genes <- cbind(res@rownames,data.frame(res@listData)) %>% 
  filter(log2FoldChange > quantile(res$log2FoldChange)[4] & padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull("res@rownames")

#remove mito genes and ribo genes
resistant_vs_control_de_genes <- resistant_vs_control_de_genes[!grepl("^MT", resistant_vs_control_de_genes)]
resistant_vs_control_de_genes <- resistant_vs_control_de_genes[!grepl("^RPS", resistant_vs_control_de_genes)]
resistant_vs_control_de_genes <- resistant_vs_control_de_genes[!grepl("^RPL", resistant_vs_control_de_genes)]

resistant_vs_control_de_signature <- resistant_vs_control_de_genes[1:signature_length]

saveRDS(resistant_vs_control_de_signature, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistant_vs_control_de_signature_", signature_length, ".rds"))

################################################################################
#Create common Resistance Signature
################################################################################

breast_de_res <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_breast_resistant_vs_control_de_res.rds")
melanoma_de_res <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_melanoma_resistant_vs_control_de_res.rds")


breast_logfc <- cbind(breast_de_res@rownames,data.frame(breast_de_res@listData)) %>% 
  filter(log2FoldChange > quantile(breast_de_res$log2FoldChange)[4] & padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select("breast_de_res@rownames","log2FoldChange")

melanoma_logfc <- cbind(melanoma_de_res@rownames,data.frame(melanoma_de_res@listData)) %>% 
  filter(log2FoldChange > quantile(melanoma_de_res$log2FoldChange)[4] & padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select("melanoma_de_res@rownames","log2FoldChange")


colnames(breast_logfc) <- c("gene","log2FC")
colnames(melanoma_logfc) <- c("gene","log2FC")

common_logfc <- merge(breast_logfc,melanoma_logfc, by="gene")

common_de_genes <- common_logfc %>% 
  mutate("avg_logFC" = rowMeans(select(.,starts_with("log")))) %>% 
  arrange(desc(avg_logFC)) %>% 
  pull(gene)

#remove mito genes
common_de_genes <- common_de_genes[!grepl("^MT", common_de_genes)]
common_de_genes <- common_de_genes[!grepl("^RPS", common_de_genes)]
common_de_genes <- common_de_genes[!grepl("^RPL", common_de_genes)]

common_de_signature <- common_de_genes[1:signature_length]

saveRDS(common_de_signature, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistant_vs_control_de_signature_", signature_length,".rds"))



raj_breast_500 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_resistant_vs_control_de_signature_500.rds")
raj_breast_200 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_resistant_vs_control_de_signature_200.rds")

raj_melanoma_500 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_melanoma_resistant_vs_control_de_signature_500.rds")
raj_melanoma_200 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_melanoma_resistant_vs_control_de_signature_200.rds")

raj_common_500 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistant_vs_control_de_signature_500.rds")
raj_common_200 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistant_vs_control_de_signature_200.rds")



raj_resistance_signatures <- list("raj_breast_500" = raj_breast_500, "raj_breast_200" = raj_breast_200,
                                 "raj_melanoma_500" = raj_melanoma_500, "raj_melanoma_200" = raj_melanoma_200,
                                 "raj_common_500" = raj_common_500, "raj_common_200" = raj_common_200)



saveRDS(raj_resistance_signatures,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_resistance_signatures.rds")

 
# refined_raj_markers <- list()
# 
# for(curr_markers in cluster_de_markers){
#   
#   curr_refined <- list(curr_markers[curr_markers %in% resistant_vs_control_de_genes])
#   
#   if(lengths(curr_refined) > 200){
#     curr_refined <- list(curr_refined[[1]][1:200])
#   }
#   
#   refined_raj_markers <- append(refined_raj_markers, curr_refined)
# }
# 
# 
# names(refined_raj_markers) <- paste0("raj_cluster",1:length(refined_raj_markers),"_refined_markers")
# 
# 
# saveRDS(refined_raj_markers, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistance_refined_cluster_markers.rds"))
# 
# refined_raj_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistance_refined_cluster_markers.rds"))
