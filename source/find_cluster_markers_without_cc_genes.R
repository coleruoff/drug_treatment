setwd("/data/ruoffcj/projects/drug_treatment/")
library(presto)
library(Seurat)
source("source/read_in_all_cell_lines.R")

cell_lines <- c("A549","K562","MCF7")

cell_cycle_genes <- unlist(cc.genes)

cluster_signatures <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- all_data[[curr_cell_line]]

  #Find global cluster markers for current cell line without cell cycle genes considered
  de_res <- FindAllMarkers(data, genes.use = setdiff(rownames(data),cell_cycle_genes))
  
  saveRDS(de_res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_cluster_de_without_cc_genes.rds"))
  # de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_cluster_de_without_cc_genes.rds"))

}


