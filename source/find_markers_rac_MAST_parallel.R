setwd("/data/ruoffcj/projects/drug_treatment/")
library(presto)
library(Seurat)
library(foreach)
library(doMC)
library(doParallel)
source("source/read_in_all_cell_lines.R")
set.seed(42)

cell_lines <- c("A549","K562","MCF7")

registerDoMC(cores=future::availableCores())

results_list <- foreach(i=cell_lines) %dopar% {
  library(Seurat)
  
  data <- all_data[[curr_cell_line]]
  cat(i,"\n")
  
  Idents(data) <- data$rac
  
  de_res <- FindMarkers(data, ident.1 = "rac", test.use = "MAST")
  
  return(de_res)
}

stopCluster()

names(results_list) <- cell_lines

saveRDS(results_list, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_global_rac_de_MAST.rds"))





