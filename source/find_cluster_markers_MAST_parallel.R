setwd("/data/ruoffcj/projects/drug_treatment/")
library(presto)
library(Seurat)
library(foreach)
library(doMC)
library(doParallel)
# source("source/read_in_all_cell_lines.R")
set.seed(42)

cell_lines <- c("A549","K562","MCF7")

curr_cell_line <- cell_lines[3]

cat(curr_cell_line,"\n")

#Read in cell line data
# data <- all_data[[curr_cell_line]]
data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))

clusters <- as.numeric(sort(unique(data$Cluster)))

registerDoMC(cores=future::availableCores())

cat("STARTING\n")
results_list <- foreach(i=1:length(clusters)) %dopar% {
  library(Seurat)
  cat(i,"\n")
  de_res <- FindMarkers(data, ident.1 = i, test.use = "MAST")
  
  return(de_res)
}

stopCluster()

names(results_list) <- clusters

saveRDS(results_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_cluster_de_MAST.rds"))





