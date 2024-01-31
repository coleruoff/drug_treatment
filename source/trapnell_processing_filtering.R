library(Seurat)
library(monocle3)
library(tidyverse)

curr_cell_line <- "A549"

#Read in trapnell data and convert it to seurat
data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raw_data/sciPlex_data/GSM4150378_sciPlex3_", curr_cell_line, "_24hrs.RDS"))
rownames(data) <- rowData(data)$gene_short_name
data <- as.Seurat(data, data=NULL)


# Select counts and create new filtered object 
counts <- data@assays$originalexp@counts

data.obj <- CreateSeuratObject(counts=counts, project=paste0(curr_cell_line,"_trapnell"), min.cells = 100, min.features = 200)

data <- data[rownames(data) %in% rownames(data.obj),colnames(data) %in% colnames(data.obj)]

data.obj@reductions <- data@reductions
data.obj@meta.data <- data@meta.data


data.obj[["percent.mt"]] <- PercentageFeatureSet(data.obj, pattern = "^MT-")

data.obj <- NormalizeData(data.obj)

data.obj <- AddMetaData(data.obj, ifelse(data.obj$dose == 0, "pre", "post"), col.name = "treatment_stage")

Idents(data.obj) <- data.obj$Cluster


saveRDS(data.obj, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))

