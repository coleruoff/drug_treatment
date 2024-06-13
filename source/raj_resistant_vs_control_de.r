args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")

library(Seurat)
library(tidyverse)
library(DESeq2)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"


if(!file.exists(paste0(dataDirectory,"de_results/"))){
  dir.create(paste0(dataDirectory,"de_results/")) 
}

################################################################################

for(raj_resistant_cancer_type in c("breast", "melanoma")){
  cat(raj_resistant_cancer_type, "\n")
  
  # Read in resistant data
  raj_resistant <- readRDS(paste0(dataDirectory,"processed_data/diverse_clonal_fates_data/raj_resistant_", raj_resistant_cancer_type, "_processed.rds"))
  
  # Psuedo-bulk resistant data
  raj_resistant_bulk <- matrix(NA, nrow=nrow(raj_resistant[["RNA"]]$counts), ncol=nlevels(raj_resistant))
  for(i in levels(raj_resistant)){
    
    raj_resistant_bulk[,as.numeric(i)+1] <- rowSums(raj_resistant[["RNA"]]$counts[,colnames(raj_resistant)[raj_resistant$seurat_clusters == i]])
    
  }
  rownames(raj_resistant_bulk) <- rownames(raj_resistant[["RNA"]]$counts)
  colnames(raj_resistant_bulk) <- paste0("raj_resistant", 1:ncol(raj_resistant_bulk))
  
  # Read in control data
  if(raj_resistant_cancer_type == "breast"){
    raj_control_cancer_type <- "MDA"
  } else if(raj_resistant_cancer_type == "melanoma"){
    raj_control_cancer_type <- "WM989"
  }
  
  raj_control <- readRDS(paste0(dataDirectory,"processed_data/memorySeq_data/", raj_control_cancer_type,"_control_counts.rds"))
  
  colnames(raj_control) <- paste0("raj_control", 1:ncol(raj_control))
  
  # Subset both datasets to only common genes
  shared_genes <- intersect(rownames(raj_resistant_bulk),rownames(raj_control))
  raj_resistant_bulk <- raj_resistant_bulk[shared_genes,]
  raj_control <- raj_control[shared_genes,]
  
  # Join two datasets and run DESeq
  all_raj_data <- cbind(raj_resistant_bulk,raj_control)
  coldata <- matrix(c(rep("resistant",ncol(raj_resistant_bulk)), rep("control", ncol(raj_control))), ncol=1)
  colnames(coldata) <- c("condition")
  rownames(coldata) <- colnames(all_raj_data)
  coldata <- data.frame(coldata)
  
  coldata$condition <- factor(coldata$condition)
  
  
  dds <- DESeqDataSetFromMatrix(countData = all_raj_data,
                                colData = coldata,
                                design = ~ condition)
  
  keep <- rowSums(counts(dds)) >= 10
  
  dds <- dds[keep,]
  
  dds$condition <- relevel(dds$condition, ref="control")
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  res <- as.matrix(res)
  
  # Save DESeq results
  saveRDS(res, paste0(dataDirectory, "de_results/raj_", raj_resistant_cancer_type,"_vs_control_deseq_result.rds"))
}


