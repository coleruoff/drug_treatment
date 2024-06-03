setwd("/data/ruoffcj/projects/drug_treatment/")
library(presto)
library(Seurat)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures"
args = commandArgs(trailingOnly=TRUE)
dataDirectory <- args[1]

################################################################################

args <- commandArgs(trailingOnly=TRUE)

curr_cell_line <- args[2]
curr_cluster <- args[3]

cat(curr_cell_line,"\n")

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

#Read in cell line data
data <- all_data[[curr_cell_line]]

# Retain genes if only expressed in > 1% of cells
genes_to_keep <- apply(data@assays$RNA$data, 1, FUN = function(x) sum(x > 0) > (.01*length(x)))
data <- data[genes_to_keep]

Idents(data) <- data$Cluster

de_res <- FindMarkers(data, ident.1 = curr_cluster, test.use = "MAST")

saveRDS(de_res, paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_",curr_cluster,"_de_MAST.rds"))

