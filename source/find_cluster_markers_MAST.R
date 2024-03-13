setwd("/data/ruoffcj/projects/drug_treatment/")
library(presto)
library(Seurat)
source("source/read_in_all_cell_lines.R")
set.seed(42)


args <- commandArgs(trailingOnly=TRUE)

curr_cell_line <- args[1]


cat(curr_cell_line,"\n")

#Read in cell line data
data <- all_data[[curr_cell_line]]

de_res <- FindAllMarkers(data, test.use = "MAST")

saveRDS(de_res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_cluster_de_MAST.rds"))







