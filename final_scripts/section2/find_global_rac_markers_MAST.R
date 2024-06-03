setwd("/data/ruoffcj/projects/drug_treatment/")
library(presto)
library(Seurat)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
args = commandArgs(trailingOnly=TRUE)
dataDirectory <- args[1]

################################################################################

args <- commandArgs(trailingOnly=TRUE)

curr_cell_line <- args[1]

cat(curr_cell_line,"\n")

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

#Read in cell line data
data <- all_data[[curr_cell_line]]

Idents(data) <- data$rac

de_res <- FindAllMarkers(data, test.use = "MAST")

saveRDS(de_res, paste0(dataDirectory, "de_results/", curr_cell_line, "_global_rac_de_MAST.rds"))







