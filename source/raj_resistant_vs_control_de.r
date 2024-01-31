library(Seurat)
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)

raj_resistant_cancer_type <- "melanoma"

raj_resistant <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/diverse_clonal_fates_data/raj_resistant_", raj_resistant_cancer_type, "_processed.rds"))

sum(raj_resistant$seurat_clusters == 1)

# rowSums(raj_resistant[["RNA"]]@data[,colnames(raj_resistant)[raj_resistant$seurat_clusters == 1]])

# str(raj_resistant)

raj_resistant_bulk <- matrix(NA, nrow=nrow(raj_resistant[["RNA"]]@data), ncol=nlevels(raj_resistant))

for(i in levels(raj_resistant)){
    
    raj_resistant_bulk[,as.numeric(i)] <- rowSums(raj_resistant[["RNA"]]@data[,colnames(raj_resistant)[raj_resistant$seurat_clusters == i]])
    
    
}

rownames(raj_resistant_bulk) <- rownames(raj_resistant[["RNA"]]@data)
colnames(raj_resistant_bulk) <- paste0("raj_resistant", 1:ncol(raj_resistant_bulk))

raj_resistant_bulk

if(raj_resistant_cancer_type == "breast"){
    raj_control_cancer_type <- "MDA"
} else if(raj_resistant_cancer_type == "melanoma"){
    raj_control_cancer_type <- "WM989"
}

raj_control <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/memorySeq_data/", raj_control_cancer_type,"_control_counts.rds"))

str(raj_control)

colnames(raj_control) <- paste0("raj_control", 1:ncol(raj_control))
dim(raj_control)

shared_genes <- intersect(rownames(raj_resistant_bulk),rownames(raj_control))

raj_resistant_bulk <- raj_resistant_bulk[shared_genes,]
raj_control <- raj_control[shared_genes,]

all_raj_data <- cbind(raj_resistant_bulk,raj_control)

all_raj_data

coldata <- matrix(c(rep("resistant",ncol(raj_resistant_bulk)), rep("control", ncol(raj_control))), ncol=1)
colnames(coldata) <- c("condition")
rownames(coldata) <- colnames(all_raj_data)
coldata <- data.frame(coldata)

coldata$condition <- factor(coldata$condition)

# coldata

dds <- DESeqDataSetFromMatrix(countData = all_raj_data,
                      colData = coldata,
                      design = ~ condition)

dds

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref="control")

dds <- DESeq(dds)

res <- results(dds)
# res

res <- as.matrix(res)
saveRDS(res, paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/raj_", raj_resistant_cancer_type,"_vs_control_deseq_result.rds"))

all_de_genes <- cbind(res@rownames,data.frame(res@listData)) %>% 
    filter(log2FoldChange > 0 & padj < 0.05) %>% 
    arrange(padj) %>% 
    pull("res@rownames")

all_de_genes <- all_de_genes[!grepl("^MT", all_de_genes)]
all_de_genes <- all_de_genes[!grepl("^RPS", all_de_genes)]
all_de_genes <- all_de_genes[!grepl("^RPL", all_de_genes)]

raj_total_resistance_geneset <- all_de_genes[1:500]

saveRDS(raj_total_resistance_geneset, paste0("/data/CDSL_hannenhalli/Cole/genesets/raj_",raj_resistant_cancer_type, "_resistant_de_signature_500.rds"))

resistant_sample_names <- colnames(raj_resistant_bulk)

resistant_cluster_genesets <- list()
for(curr_resistant_sample in resistant_sample_names){
    curr_geneset <- names(sort(raj_resistant_bulk[all_de_genes,curr_resistant_sample], decreasing = T)[1:100])
    
    resistant_cluster_genesets <- append(resistant_cluster_genesets, list(curr_geneset))
}

resistant_cluster_genesets <- append(resistant_cluster_genesets,list(raj_total_resistance_geneset))

names(resistant_cluster_genesets) <- c(paste0(resistant_sample_names, "_geneset"), "raj_total_resistance_geneset")

resistant_cluster_genesets

saveRDS(resistant_cluster_genesets, paste0("/data/CDSL_hannenhalli/Cole/genesets/raj_",raj_resistant_cancer_type, "_resistant_de_genesets.rds"))



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
    
    raj_cluster_de_genesets <- append(raj_cluster_de_genesets, list(all_de_genes[1:200]))
}
                                      
names(raj_cluster_de_genesets) <- paste0(colnames(raj_resistant_bulk), "_de_genes")

saveRDS(raj_cluster_de_genesets, paste0("/data/CDSL_hannenhalli/Cole/genesets/raj_",raj_resistant_cancer_type,"_cluster_de_genesets_200.rds"))
