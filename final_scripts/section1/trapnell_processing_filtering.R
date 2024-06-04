args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

library(Seurat)
library(monocle3)
library(tidyverse)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

if(!file.exists(paste0(dataDirectory,"processed_data/sciPlex_data/"))){
 dir.create(paste0(dataDirectory,"processed_data/sciPlex_data/")) 
}

cell_lines <- c("A549","K562","MCF7")

for(curr_cell_line in cell_lines){
  
  #Read in Trapnell data and convert it to Seurat object
  data <- readRDS(paste0(dataDirectory, "raw_data/sciPlex_data/GSM4150378_sciPlex3_", curr_cell_line, "_24hrs.RDS"))
  rownames(data) <- rowData(data)$gene_short_name
  data <- as.Seurat(data, data=NULL)
  
  # Select counts and create new filtered object 
  counts <- data@assays$originalexp@counts
  
  data.obj <- CreateSeuratObject(counts=counts, project=paste0(curr_cell_line,"_trapnell"), min.cells = 100, min.features = 200)
  
  data <- data[rownames(data) %in% rownames(data.obj),colnames(data) %in% colnames(data.obj)]
  
  # data.obj@assays$RNA$counts
  data.obj@reductions <- data@reductions
  data.obj@meta.data <- data@meta.data
  colnames(data.obj@meta.data) <- sapply(colnames(data.obj@meta.data), FUN = function(x) gsub("_originalexp","_RNA", x))
  
  data.obj[["percent.mt"]] <- PercentageFeatureSet(data.obj, pattern = "^MT-")
  
  data.obj <- NormalizeData(data.obj)
  
  data.obj <- AddMetaData(data.obj, ifelse(data.obj$dose == 0, "pre", "post"), col.name = "treatment_stage")
  
  Idents(data.obj) <- data.obj$Cluster
  
  saveRDS(data.obj, paste0(dataDirectory,"processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
}